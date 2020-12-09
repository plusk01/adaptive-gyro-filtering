/**
 * @file adaptnotch.h
 * @brief Adaptive notch filtering algorithm
 * @author Parker Lusk <parkerclusk@gmail.com>
 * @date 8 Dec 2020
 */

#pragma once

#include <algorithm>
#include <complex>
#include <memory>
#include <vector>

#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>

namespace adaptnotch {

  /**
   * @brief      Ring buffer implementation for applying FFT to chunks of data
   */
  class RingBuffer
  {
  public:
    /**
     * @brief      RingBuffer constructor. Values initialize to zero.
     *
     * @param[in]  N     Fixed size of buffer
     */
    RingBuffer(size_t N) { data_ = Eigen::VectorXd::Zero(N); }
    ~RingBuffer() = default;

    /**
     * @brief      Push a new element onto the buffer.
     *
     * @param[in]  x     Element to add
     */
    void add(double x)
    {
      data_(idx)= x;
      idx = (idx + 1) % data_.size();
    }

    /**
     * @brief      Copies the internal buffer into a sequential buffer.
     *
     * @return     Vector of elements ordered from last (0) to first (N).
     */
    Eigen::VectorXd sequentialView() const
    {
      Eigen::VectorXd v = Eigen::VectorXd::Zero(data_.size());
      const size_t m = data_.size() - idx;
      // push the first m elements onto v
      for (size_t i=0; i<m; ++i) {
        v(i) = data_(idx+i);
      }
      // push the remaining (wrapped) elements onto end of v
      if (idx > 0) {
        for (size_t i=m; i<data_.size(); ++i) {
          v(i) = data_(i-m);
        }
      }
      return v;
    }

    /**
     * @brief      View internal ring buffer
     *
     * @return     View of internal ring buffer
     */
    const Eigen::VectorXd& data() const { return data_; }

  private:
    size_t idx = 0; ///< index of next spot to add element at
    Eigen::VectorXd data_; ///< internal ring buffer
  };

  /**
   * @brief      Digital biquadratic notch filter implementation. Uses
   *             Direct Form I (DF1) implementation which allows updating
   *             the coefficients of the filter on the fly.
   */
  class BiquadNotch
  {
  public:
    /**
     * @brief      BiquadNotch constructor.
     *
     * @param[in]  fc    Center frequency of notch
     * @param[in]  fs    Sample rate of data filter will be applied to
     * @param[in]  Q     Filter quality factor (higher == tighter bandwidth)
     */
    BiquadNotch(double fc, double fs, double Q) { init(fc, fs, Q); }
    ~BiquadNotch() = default;

    /**
     * @brief      Updates the filter with new specifications
     *
     * @param[in]  fc    Center frequency of notch
     * @param[in]  fs    Sample rate of data filter will be applied to
     * @param[in]  Q     Filter quality factor (higher == tighter bandwidth)
     */
    void update(double fc, double fs, double Q)
    {
      // backup state
      double x1 = x1_, x2 = x2_;
      double y1 = y1_, y2 = y2_;

      // recreate filter with new specs
      init(fc, fs, Q);

      // restore state
      x1_ = x1; x2_ = x2;
      y1_ = y1; y2_ = y2;
    }

    /**
     * @brief      Applies the filter to the input sample. Uses a DF1
     *             implementation to support dynamically changing filter specs.
     *
     * @param[in]  x     Input sample
     *
     * @return     Filtered output sample
     */
    double apply(double x)
    {
      // apply filter using direct form 1 (DF1)
      const double y = b0_ * x  +  b1_ * x1_  +  b2_ * x2_
                     -  (a1_ * y1_  +  a2_ * y2_);

      // shift feedback delay lines
      x2_ = x1_;
      x1_ = x;

      // shift feedforward delay lines
      y2_ = y1_;
      y1_ = y;

      return y;
    }

  private:
    // \brief Filter parameters
    double fc_, fs_, Q_;  ///< center freq, sample freq, quality factor
    double b0_, b1_, b2_; ///< num coeffs
    double a1_, a2_;      ///< den coeffs

    // \brief Filter state
    double x1_, x2_; ///< feedback delay elements
    double y1_, y2_; ///< feedforward delay elements

    /**
     * @brief      Create a biquad notch filter
     *
     * @param[in]  fc    Center frequency of notch
     * @param[in]  fs    Sample rate of data filter will be applied to
     * @param[in]  Q     Filter quality factor (higher == tighter bandwidth)
     */
    void init(double fc, double fs, double Q)
    {
      // normalized frequency in [0, pi]
      static constexpr double TWO_PI = 2 * 3.14159265358979323846;
      const double omega = TWO_PI * (fc / fs);

      const double sn = std::sin(omega);
      const double cs = std::cos(omega);
      const double alpha = sn / (2 * Q);

      // keep around just for fun
      fc_ = fc;
      fs_ = fs;
      Q_ = Q;

      // notch biquad setup
      const double b0 = 1;
      const double b1 = -2 * cs;
      const double b2 = 1;
      const double a0 = 1 + alpha;
      const double a1 = -2 * cs;
      const double a2 = 1 - alpha;

      // normalize into standard biquad form (a0 == 1)
      b0_ = b0 / a0;
      b1_ = b1 / a0;
      b2_ = b2 / a0;
      a1_ = a1 / a0;
      a2_ = a2 / a0;

      // initialize delay elements
      x1_ = x2_ = 0;
      y1_ = y2_ = 0;
    }
  };

  class AdaptiveNotch
  {
  public:
    /**
     * @brief      Adaptive notch algorithm parameters 
     */
    struct Params
    {
      int Fs = 500; ///< sample rate of input data
      int NFFT = 128; ///< buffer and FFT length
      int dual_notch_width_percent = 8; ///< 0 width is single notch
      int Q = 360; ///< bandwidth of notch
      int min_hz = 60; ///< lower bound of peak detection
      int max_hz = 200; ///< upper bound of peak detection
      int start_bin = 2; ///< skip bins before this in peak search (e.g., DC)
    };

  public:
    AdaptiveNotch(const Params& params);
    ~AdaptiveNotch() = default;
    
    double apply(double x);

  private:
    Params params_; ///< instance parameters
    Eigen::FFT<double> fft_; ///< fft object
    std::unique_ptr<BiquadNotch> notch1_, notch2_; ///< dual notch filters

    // \brief Parameter initialization
    int fft_bin_count_; ///< num useful bins for real input, excl. Nyquist bin
    double notch1_ctr_, notch2_ctr_; ///< dual notch scale factors
    double Q_; ///< Filter quality factor
    int min_hz_, max_hz_; ///< search range for peaks
    int fftFs_; ///< sample rate of downsampled input data
    int max_samples_; ///< number of samples to accumulate before downsampling
    double fres_; ///< frequency resolution
    int start_bin_; ///< start peak search at this bin
    Eigen::VectorXd window_; ///< window for tapering data for FFT

    // \brief State
    double peakFreq_; ///< estimated frequency of noise peak in specified range

    double input_accumulator_; ///< accumulator for downsampling to fft Fs
    double input_samples_; ///< num samples in accumulator

    RingBuffer buffer_;

    void reset();
    Eigen::VectorXd windowHann(int N);
    double findFreqPeak();

  };

  namespace utils {
    template<typename T>
    T clamp(T v, T lb, T ub) {
      if (v < lb) v = lb;
      if (v > ub) v = ub;
      return v;
    }
  } // ns utils

} // ns adaptnotch