// Copyright 2003-2005 Scott Jackson
// Copyright 2008-2012 Muhammad S. A. Zilany
// Copyright 2019 Google LLC
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file or at
// https://developers.google.com/open-source/licenses/bsd

// Fractional Gaussian noise and Brownian motion generators. Brownian motion is
// an important continuous-time stochastic process that serves as a
// continuous-time analog to the simple symmetric random walk on the one hand,
// and shares fundamental properties with the Poisson counting process on the
// other hand.
//
// The APIs in this library are mostly C++ reimplementations of the Matlab
// originals. Please see fractional_gaussian_noise.README file in this directory
// for pointers to the original code, documentation and licenses.
//
// See:
//   - https://en.wikipedia.org/wiki/Fractional_Brownian_motion
//   - Mandelbrot, B. B., & Van Ness, J. W. (1968)., "Fractional Brownian
//     motions, fractional noises and applications", SIAM review, 10(4),
//     422-437.
//   - http://www.columbia.edu/~ks20/FE-Notes/4700-07-Notes-BM.pdf

#ifndef THIRD_PARTY_AUDITION_MODELS_BRUCE_CARNEY_FRACTIONAL_GAUSSIAN_NOISE_H_
#define THIRD_PARTY_AUDITION_MODELS_BRUCE_CARNEY_FRACTIONAL_GAUSSIAN_NOISE_H_

#include "Eigen/Core"

namespace eidos {
namespace audition {

// FFGN: A fast (exact) fractional Gaussian noise and Brownian motion generator.
// Returns a vector containing a sequence of fractional Gaussian noise or
// fractional Brownian motion. The generation process uses an FFT which makes
// it very fast.
//
// The length of the output sequence is given by <output_length>. Hurst index
// of the resultant noise (<hurst_index>) must be in the range (0 < H <= 2). For
// 0 < H <= 1, the output will be fractional Gaussian noise with Hurst index H.
// For 1 < H <= 2, the output will be fractional Brownian motion with Hurst index
// H-1. Either way, the power spectral density of the output will be nominally
// proportional to 1/f^(2H-1). The mean and standard deviation of the resulting
// noise are specified by <mu> and <sigma> parameters, respectively. A sensible
// choice of these parameters is zero mean and unit standard deviation.
//
// This method is based on an embedding of the covariance matrix in a circulant
// matrix.
//
// References:
//   - Davies, R.B. and Harte, D.S. (1987): "Tests for Hurst effect.",
//     Biometrika, 74(1), pp. 95-101.
//   - Beran, J. (2017): "Statistics for long-memory processes.", Routledge.
//   - Bardet, J.M. and Bertrand, P. (2002): "Detecting abrupt change on the
//     Hurst parameter of a multi-scale fractional Brownian motion with
//     applications.", Internation Meeting on Self-Similarity and Applications,
//     pp. 27-30.
//
// Original implementation: Copyright (c) 2003-2005 by B. Scott Jackson.
Eigen::VectorXd JacksonFractionalGaussianNoise(int output_length,
                                               double hurst_index,
                                               double mu,
                                               double sigma);

}  // namespace audition
}  // namespace eidos

#endif  // THIRD_PARTY_AUDITION_MODELS_BRUCE_CARNEY_FRACTIONAL_GAUSSIAN_NOISE_H_

// Local Variables:
// mode: c++
// End:
