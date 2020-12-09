/**
 * @file adaptnotch.h
 * @brief Adaptive notch filtering algorithm
 * @author Parker Lusk <parkerclusk@gmail.com>
 * @date 8 Dec 2020
 */

#pragma once

#include <complex>
#include <vector>

#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>

namespace adaptnotch {

  class RingBuffer
  {
  public:
    RingBuffer(size_t N) { data_ = Eigen::VectorXd::Zero(N); }
    ~RingBuffer() = default;

    void add(double x)
    {
      data_(idx)= x;
      idx = (idx + 1) % data_.size();
    }

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

    const Eigen::VectorXd& data() const { return data_; }

  private:
    size_t idx = 0; ///< index of next spot to add element at
    Eigen::VectorXd data_;
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

    // \brief Parameter initialization
    int fft_bin_count_; ///< num useful bins for real input, excl. Nyquist bin
    double notch1_ctr_, notch2_ctr_; ///< dual notch scale factors
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