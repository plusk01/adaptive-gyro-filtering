/**
 * @file main.cpp
 * @brief Entry-point for dynamic notch filter example
 * @author Parker Lusk <parkerclusk@gmail.com>
 * @date 8 Dec 2020
 */

#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>

#include <plot.hpp>

#include <adaptnotch/adaptnotch.h>

#include "csv.h"

static constexpr int DATUMS = 4; // number of columns to be extracted from CSV
using Data = Eigen::Matrix<double, Eigen::Dynamic, DATUMS>;

void usage(int argc, char const *argv[])
{
  std::cout << argv[0] << " <input csv data> [axis] [plot]" << std::endl << std::endl;
  std::cout << "\tRun adaptive notching algorithm on gyro data stored in CSV.";
  std::cout << std::endl << "\t";
  std::cout << "CSV file expected to have been generated from sfpro or to be ";
  std::cout << std::endl << "\t";
  std::cout << "in same format. A gyro axis to analyze may be specified as";
  std::cout << std::endl << "\t";
  std::cout << "(1, 2, 3) which corresponds to axes (x, y, z).";
  std::cout << std::endl << std::endl;
  std::cout << "\tIf a 3rd argument is specified (e.g., '1'), an ASCII plot";
  std::cout << std::endl << "\t";
  std::cout << "of the pre- and post-spectrum is shown in the terminal.";
  std::cout << std::endl << std::endl;
}

// ----------------------------------------------------------------------------

Data parseCSV(const std::string& file)
{

  // count number of entries (estimate)
  std::ifstream ifile(file);
  const int N = std::count(std::istreambuf_iterator<char>(ifile),
                            std::istreambuf_iterator<char>(), '\n');
  ifile.close();

  io::CSVReader<DATUMS> in(file);
  in.next_line(); // ignore "dsp clock" message

  // we only care about these four (DATUMS) columns
  in.read_header(io::ignore_extra_column,
                  "timestamp(us)", "ang_x", "ang_y", "ang_z");

  Data D = Data::Zero(N, DATUMS);
  int i = 0;
  int time_us;
  double wx, wy, wz;
  while (in.read_row(time_us, wx, wy, wz)) {
    D.row(i++) << time_us*1e-6, wx, wy, wz;
  }

  // resize to however many valid entries there were
  D.conservativeResize(i, DATUMS);

  return D;
}

// ----------------------------------------------------------------------------

void plotResults(const adaptnotch::AdaptiveNotch& filter,
            const Eigen::VectorXd& gyro, const Eigen::VectorXd& gyrof, int n)
{
  static const int N = filter.params().NFFT;
  static Eigen::FFT<double> fft;
  static plot::TerminalInfo term;
  static bool init = false;
  if (!init) {
    term.detect();
    fft.SetFlag(Eigen::FFT<double>::Unscaled);
    fft.SetFlag(Eigen::FFT<double>::HalfSpectrum);
    init = true;
  }
  constexpr float ymax = 0.2f; // arbitrary FFT mag scaling
  static plot::RealCanvas<plot::BrailleCanvas> prefilter({ { 0.0f, ymax }, { N/2.0f, 0.0f } }, plot::Size(60, 10), term);
  static plot::RealCanvas<plot::BrailleCanvas> postfilter({ { 0.0f, ymax }, { N/2.0f, 0.0f } }, plot::Size(60, 10), term);

  // Build block layout
  auto layout =
      plot::alignment(
          { term.size().x, 0 },
          plot::margin(
                  plot::vbox(
                      plot::frame(u8"Original Spectrum", plot::Align::Center, &prefilter),
                      plot::frame(u8"Filtered Spectrum", plot::Align::Center, &postfilter))));

  // pre-filter spectrum
  const Eigen::VectorXd Y = filter.spectrum() / N;

  // select N most recent filtered measurements
  const size_t s = (n-N<0) ? 0 : n-N;
  const Eigen::VectorXd yf = gyrof.segment(s,N);

  // post-filter spectrum
  Eigen::VectorXcd Yfc;
  fft.fwd(Yfc, yf);
  const Eigen::VectorXd Yf = Yfc.array().abs() / N;

  // Plot spectrum pre-filtering
  prefilter.clear();
  for (size_t i=1; i<N/2; i++) {
    prefilter.path(plot::palette::royalblue,{{static_cast<float>(i-1), static_cast<float>(Y(i-1))},
                                          {static_cast<float>(i), static_cast<float>(Y(i))}});
  }

  // Plot spectrum post-filtering
  postfilter.clear();
  for (size_t i=1; i<N/2; i++) {
    postfilter.path(plot::palette::royalblue,{{static_cast<float>(i-1), static_cast<float>(Yf(i-1))},
                                          {static_cast<float>(i), static_cast<float>(Yf(i))}});
  }

  for (auto const& line: layout)
    std::cout << term.clear_line() << line << std::endl;
  std::cout << term.move_up(layout.size().y) << std::flush;
}

// ----------------------------------------------------------------------------

int main(int argc, char const *argv[])
{

  int axis = 1;     ///< x, y, or z axis of gyro to analyze
  std::string file; ///< input data from IMU
  bool shouldPlot = false; ///< show FFT plots in terminal

  if (argc < 2) {
    std::cerr << "Not enough input arguments." << std::endl << std::endl;
    usage(argc, argv);
    return -1;
  } else if (argc >= 2) {
    file = std::string(argv[1]);
  }

  if (argc >= 3) {
    axis = std::stoi(argv[2]);
    if (axis < 1 || axis > 3) axis = 1;
  }

  if (argc >= 4) shouldPlot = true;

  //
  // Process raw IMU data
  //

  Data D = parseCSV(file);

  const Eigen::VectorXd diff = D.col(0).bottomRows(D.rows()-1) - D.col(0).topRows(D.rows()-1);
  const double Ts = diff.mean();
  const double Fs = 1./Ts;
  const int N = D.rows();

  //
  // Adaptive notch filter setup
  //

  adaptnotch::AdaptiveNotch::Params params;
  adaptnotch::AdaptiveNotch filter(params);

  //
  // Main loop - simulated gyro sampling
  //

  Eigen::VectorXd gyrof = Eigen::VectorXd::Zero(N);
  const auto start = std::chrono::steady_clock::now();


  for (size_t n=0; n<N; n++) {
    const double gyro = D(n, axis);
    gyrof(n) = filter.apply(gyro);
    if (shouldPlot) plotResults(filter, D.col(axis), gyrof, n);
  }


  //
  // Timing stats
  //

  const auto end = std::chrono::steady_clock::now();
  const double duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() * 1e-6;

  const double timu = D(N-1,0) - D(0,0);

  std::cout << "Processed " << timu << " seconds (" << N << " samples) of IMU";
  std::cout << " data in " << duration << " seconds" << std::endl;
  std::cout << "Real-time factor: " << timu / duration << std::endl;
  std::cout << "Estimated peak freq: " << filter.peakFreq() << std::endl;

  //
  // Write data to file
  //

  Eigen::MatrixXd out(N, 3);
  out << Eigen::VectorXd::LinSpaced(N, 0, N*Ts), D.col(axis), gyrof;

  std::ofstream of("data_processed.txt");
  of << out << std::endl;
  of.close();

  return 0;
}
