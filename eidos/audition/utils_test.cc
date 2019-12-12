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

// Sanity checks for utils library.

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

#include "eidos/audition/utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

using testing::AllOf;
using testing::Each;
using testing::Ge;
using testing::Le;

namespace eidos {
namespace audition {

TEST(UtilsTest, CheckRandomUniformVector) {
  constexpr int kMaxLengthTrials = 10;
  constexpr int kMaxNumRanges = 100;
  for (int l = 1; l <= kMaxLengthTrials; ++l) {
    const int n = l * 100;
    for (int i = 1; i <= kMaxNumRanges; ++i) {
      const std::vector<double> vec = GetRandomUniformVector(
          n, -i, i);
      EXPECT_EQ(n, vec.size());
      EXPECT_THAT(vec, Each(AllOf(Ge(-i), Le(i))));
    }
  }
}

TEST(UtilsTest, CheckRandomNormalVector) {
  constexpr int kVecLen = 100000;
  constexpr double kEpsilon = 2E-2;
  const std::vector<double> vec = GetRandomNormalVector(
      kVecLen, 0.0 /* zero mean */, 1.0 /* unit standard deviation */);
  EXPECT_EQ(kVecLen, vec.size());
  const double mean = std::accumulate(vec.begin(), vec.end(), 0.0) / vec.size();
  EXPECT_NEAR(0.0, std::fabs(mean), kEpsilon);
  const double sum_squares = std::inner_product(
      vec.begin(), vec.end(), vec.begin(), 0.0,
      [](double x, double y) { return x + y; },
      [mean](double x, double y) { return (x - mean) * (y - mean); });
  const double stddev = sum_squares / (vec.size() - 1);
  EXPECT_NEAR(1.0, stddev, kEpsilon);
}

}  // namespace audition
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
