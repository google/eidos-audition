// Copyright 2019 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    https://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// Checks that fractional Gaussian and Brownian motion APIs behave.

#include <algorithm>
#include <cmath>
#include <numeric>

#include "Eigen/Core"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "third_party/audition/models/bruce_carney/fractional_gaussian_noise.h"

using testing::Each;
using testing::Ne;

namespace eidos {
namespace audition {

TEST(FractionalGaussianNoiseTest, JacksonWhiteGaussianNoise) {
  // Check a case of pure white Gaussian noise.
  constexpr int kOutputLength = 100000;
  constexpr double kEpsilon = 2E-1;
  Eigen::VectorXd output = JacksonFractionalGaussianNoise(
      kOutputLength, 0.5 /* Hurst index */,
      0.0 /* mean */, 1.0 /* standard deviation */);
  EXPECT_EQ(kOutputLength, output.size());
  std::vector<double> vec(output.data(), output.data() + output.size());
  EXPECT_THAT(vec, Each(Ne(0.0)));
  double mean = std::accumulate(vec.begin(), vec.end(), 0.0) / vec.size();
  EXPECT_NEAR(0.0, mean, kEpsilon);

  // Same as above, but Brownian motion corresponding to white Gaussian
  // noise. Compute cumulative subtraction to get at the white noise.
  output = JacksonFractionalGaussianNoise(
      kOutputLength, 1.5 /* Hurst index */,
      0.0 /* mean */, 1.0 /* standard deviation */);
  EXPECT_EQ(kOutputLength, output.size());
  double cum_sum = 0.0;
  for (int i = 0; i < kOutputLength; ++i) {
    vec[i] = output[i] - cum_sum;
    cum_sum += vec[i];
  }
  mean = std::accumulate(vec.begin(), vec.end(), 0.0) / vec.size();
  EXPECT_NEAR(0.0, mean, kEpsilon);
}

TEST(FractionalGaussianNoiseTest, JacksonFractionalGaussianNoiseOrBrownian) {
  // Generate fractional Gaussian noise.
  constexpr int kOutputLength = 30000;
  constexpr double kEpsilon = 1.5;
  for (double hurst = 0.01; hurst <= 2.0; hurst += 0.05) {
    const Eigen::VectorXd output = JacksonFractionalGaussianNoise(
        kOutputLength, hurst /* Hurst index */,
        0.0 /* mean */, 1.0 /* standard deviation */);
    EXPECT_EQ(kOutputLength, output.size());
    const std::vector<double> vec(output.data(),
                                  output.data() + output.size());
    if (hurst < 1.0) {
      // Fractional Gaussian noise. The mean should be roughly around zero and
      // the standard deviation around unity, which, given the small sample
      // size, is going to just about hold true given a large epsilon.
      const double mean = std::accumulate(vec.begin(), vec.end(), 0.0) /
          vec.size();
      EXPECT_NEAR(0.0, std::fabs(mean), kEpsilon);
      const double sum_squares = std::inner_product(
          vec.begin(), vec.end(), vec.begin(), 0.0,
          [](double x, double y) { return x + y; },
          [mean](double x, double y) { return (x - mean) * (y - mean); });
      const double stddev = sum_squares / (vec.size() - 1);
      EXPECT_NEAR(1.0, std::fabs(stddev), kEpsilon);
    }
  }
}

}  // namespace audition
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
