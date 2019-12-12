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

// Check that various psychoacoustic scale conversions are sane.

#include "eidos/audition/psychoacoustic_scale_utils.h"

#include <algorithm>
#include <cstdlib>
#include <vector>

#include "gtest/gtest.h"

namespace eidos {
namespace audition {

TEST(PsychoacousticScaleUtilsTest, CheckErb) {
  // Checks unit conversions.
  constexpr double kSourceHz = 2000;  // 2 kHz
  constexpr double kTargetErb = 240.578;
  constexpr double kEpsilon = 1E-5;
  EXPECT_NEAR(kTargetErb, erb::HzToErb(kSourceHz), kEpsilon);
  EXPECT_NEAR(kSourceHz, erb::ErbToHz(kTargetErb), kEpsilon);
}

TEST(PsychoacousticScaleUtilsTest, CheckErbRates) {
  // Checks rate conversions. Note that the epsilon for this test is higher
  // because of the use of logarithmic and exponential functions during the
  // conversions.
  constexpr double kSourceHz = 2000;  // 2 kHz
  constexpr double kTargetErbRate = 21.1552;
  constexpr double kEpsilon = 1E-2;
  EXPECT_NEAR(kTargetErbRate, erb::HzToErbRate(kSourceHz), kEpsilon);
  EXPECT_NEAR(kSourceHz, erb::ErbRateToHz(kTargetErbRate), kEpsilon);
}

TEST(PsychoacousticScaleUtilsTest, CheckErbRange) {
  constexpr double kMinFrequency = 50;  // Hz
  constexpr double kMaxFrequency = 10000;  // Hz
  constexpr int kNumChannels = 500;
  const std::vector<double> freqs = erb::GetFrequenciesHz(
      kNumChannels, kMinFrequency, kMaxFrequency);
  // Following checks for weak ordering.
  EXPECT_TRUE(std::is_sorted(freqs.begin(), freqs.end()));
  // Following checks for strict ordering.
  for (std::size_t i = 0; i < freqs.size(); ++i) {
    if (i) {
      EXPECT_LT(freqs[i - 1], freqs[i]);
    }
  }
  EXPECT_LT(kMinFrequency, freqs[0]);
  EXPECT_GT(kMaxFrequency, freqs.back());
}

}  // namespace audition
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
