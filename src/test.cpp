
#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>
#include <complex>
#include <cmath>
#include <iostream>
#include <fstream>

unsigned const N = 1000;  //
double const Fs  = 32;    // [Hz]
double const Ts  = 1./Fs; // [s] 
const double f0  = 5;     // [Hz]

std::complex<double> f(std::complex<double> const & t) { return std::sin(2*M_PI*f0*t); }
double f(double const & t) { return std::sin(2*M_PI*f0*t); }

int main(){
    std::ofstream xrec("xrec.txt");
    Eigen::VectorXd time(N);
    Eigen::VectorXd f_values(N);
    Eigen::VectorXd freq(N);
    for(int u = 0; u < N; ++u){
        time(u) = u * Ts;
        f_values(u) = f(time(u));
        freq(u) = Fs * u / double(N);
    }

    Eigen::FFT<double> fft;
    fft.SetFlag(Eigen::FFT<double>::Unscaled);
    fft.SetFlag(Eigen::FFT<double>::HalfSpectrum);
    Eigen::VectorXcd f_freq(N);
    fft.fwd(f_freq, f_values);

    for(int u = 0; u < N; ++u){
        xrec << freq(u) << " " << std::abs(f_freq(u)) << "\n"; 
    }
};