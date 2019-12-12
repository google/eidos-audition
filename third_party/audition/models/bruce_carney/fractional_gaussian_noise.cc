// Copyright 2003-2005 Scott Jackson
// Copyright 2008-2012 Muhammad S. A. Zilany
// Copyright 2019 Google LLC
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file or at
// https://developers.google.com/open-source/licenses/bsd

#include "third_party/audition/models/bruce_carney/fractional_gaussian_noise.h"

#include <cmath>
#include <complex>
#include <vector>

#include "eidos/audition/utils.h"
#include "eidos/stubs/logging.h"
#include "unsupported/Eigen/FFT"

namespace eidos {
namespace audition {
namespace {

// Computes cumulative sum of the <input> mimicking Matlab cumsum API.
Eigen::VectorXd CumulativeSum(const Eigen::VectorXd &input) {
  const int num_elements = input.size();
  Eigen::VectorXd output(num_elements);
  double cum_sum = 0.0;
  for (int i = 0; i < num_elements; ++i) {
    cum_sum += input[i];
    output[i] = cum_sum;
  }
  return output;
}

// Computes random vectors drawn from normal distribution.
Eigen::VectorXd RandomNormalVector(int length) {
  const std::vector<double> vec = GetRandomNormalVector(
      length, 0.0 /* mean */, 1.0 /* standard deviation */);
  return Eigen::VectorXd::Map(vec.data(), vec.size());
}

}  // namespace

Eigen::VectorXd JacksonFractionalGaussianNoise(int output_length,
                                               double hurst_index,
                                               double mu,
                                               double sigma) {
  if (hurst_index <= 0.0 || hurst_index > 2.0) {
    GOOGLE_LOG(FATAL) << "The Hurst index must be in the interval (0, 2]";
  }
  if (sigma <= 0.0) {
    GOOGLE_LOG(FATAL) << "Standard deviation must be greater than zero";
  }

  // Determine whether Gaussian noise or Brownian motion should be produced.
  bool brownian_motion = false;
  if (hurst_index > 1.0) {
    brownian_motion = true;
    hurst_index -= 1.0;
  }

  // Generate the outputs.
  Eigen::VectorXd output;
  if (hurst_index == 0.5) {
    // Gaussian noise is equivalent to white Gaussian noise.
    output = RandomNormalVector(output_length);
  } else {
    const double n_fft = std::pow(2.0, std::ceil(
        std::log2(2.0 * (output_length - 1))));
    const int n_fft_half = std::round(n_fft / 2.0);
    Eigen::ArrayXd k(2 * n_fft_half);
    for (int i = 0; i <= n_fft_half; ++i) {
      k[i] = i;
    }
    for (int i = 1; i < n_fft_half; ++i) {
      k[n_fft_half + i] = n_fft_half - i;
    }
    // Zmag = 0.5 .* ((k+1).^(2.*H) - 2.*k.^(2.*H) + (abs(k-1)).^(2.*H)).
    const double double_hurst = 2.0 * hurst_index;
    Eigen::ArrayXd z_mag = 0.5 * ((k + 1.0).pow(double_hurst) -
                                  2.0 * k.pow(double_hurst) +
                                  (k - 1.0).abs().pow(double_hurst));
    // Zmag = real(fft(Zmag)).
    Eigen::FFT<double> fft;
    Eigen::VectorXcd fft_out;
    fft.fwd(fft_out, z_mag.matrix());
    z_mag = fft_out.real();
    if (z_mag.minCoeff() < 0.0) {
      GOOGLE_LOG(FATAL) << "The fast Fourier transform of the circulant "
                        << "covariance had negative values.";
    }
    z_mag = z_mag.sqrt();

    // Z = Zmag.*(randn(1,Nfft) + i.*randn(1,Nfft)).
    const Eigen::VectorXd rand_r = RandomNormalVector(n_fft);
    const Eigen::VectorXd rand_i = RandomNormalVector(n_fft);
    Eigen::ArrayXcd Z(rand_r.size());;
    Z.real() = rand_r;
    Z.imag() = rand_i;
    Z *= z_mag;

    // y = real(ifft(Z)) .* sqrt(Nfft).
    Eigen::VectorXcd Z_inv;
    fft.inv(Z_inv, Z.matrix());
    output = (Z_inv.real() * std::sqrt(n_fft)).head(output_length);
  }

  // Convert to Brownian motion, if necessary.
  if (brownian_motion) {
    output = CumulativeSum(output);
  }

  // Translate by the mean and scale by standard deviation.
  return (output.array() + mu) * sigma;
}

}  // namespace audition
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
