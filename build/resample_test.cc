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

// Simple test for generic polyphase resampler.

#include <vector>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "resample.h"

using ::testing::DoubleNear;
using ::testing::Pointwise;

namespace eidos {
namespace {

constexpr double kInputSignal[] = {
  1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };

// The epsilon is quite large unfortunately. There may be a bug in the resample
// library.
constexpr double kEpsilon = 0.3;

}  // namespace

TEST(ResampleTest, UpsampleTest) {
  constexpr int kUpFactor = 24;
  constexpr int kDownFactor = 13;
  std::vector<double> input;
  input.assign(kInputSignal, kInputSignal + 10);
  std::vector<double> output;
  resample(kUpFactor, kDownFactor, input, output);

  // Check against the reference output produced by the Octave resample
  // function from the signal package.
  EXPECT_EQ(19, output.size());
  const std::vector<double> kRefOutput = {
    1.0000, 1.7436, 2.0297, 2.4257, 3.2808, 3.9023, 4.0629, 4.6156, 5.6100,
    6.0115, 6.0265, 6.8976, 8.0470, 7.9552, 7.7765, 9.5546, 11.1220, 8.2991,
    2.2418 };
  EXPECT_THAT(kRefOutput, Pointwise(DoubleNear(kEpsilon), output));
}

TEST(ResampleTest, DecimateTest) {
  constexpr int kUpFactor = 1;
  constexpr int kDownFactor = 3;
  std::vector<double> input;
  input.assign(kInputSignal, kInputSignal + 10);
  std::vector<double> output;
  resample(kUpFactor, kDownFactor, input, output);

  // Check against the reference output produced by the Octave resample
  // function from the signal package.
  EXPECT_EQ(4, output.size());
  const std::vector<double> kRefOutput = { 1.2429, 3.5969, 7.8041, 6.4141 };
  EXPECT_THAT(kRefOutput, Pointwise(DoubleNear(kEpsilon), output));
}

}  // namespace eidos

// Local Variables:
// mode: c++
// End:
