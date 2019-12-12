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

// A very simple test for the Zhang and Carney's spike generator.

#include <vector>

#include "Eigen/Core"
#include "gtest/gtest.h"
#include "third_party/audition/models/bruce_carney/zhang2001_spike_generator.h"

namespace eidos {
namespace audition {

TEST(Zhang2001SpikeGeneratorTest, MultiFiberRandomStimulus) {
  constexpr int kNumChannels = 100;
  constexpr int kNumSamples = 5000;  // Number of samples to test.
  constexpr double kSampleRate = 16000;  // 16 kHz
  constexpr double kMaxRepeats = 5;  // Maximum number of trials.
  const std::vector<int> kMinNumSpikes = {
    500, 1000, 1500, 2000, 3000 };

  const Eigen::ArrayXXd input =
      Eigen::ArrayXXd::Random(kNumChannels, kNumSamples) * 100.0;
  for (int n = 1; n <= kMaxRepeats; ++n) {
    Eigen::ArrayXXd spikes;
    Zhang2001SpikeGenerator(input,
                            1.0 / kSampleRate,  // Sample period (sec).
                            n, &spikes);
    EXPECT_EQ(kNumChannels, spikes.rows());
    EXPECT_EQ(kNumSamples, spikes.cols());
    EXPECT_FALSE(spikes.isZero());
    EXPECT_EQ(0.0, spikes.minCoeff());
    EXPECT_EQ(1.0, spikes.maxCoeff());

    // Expect these lower bounds on the number of spikes.
    EXPECT_LT(kMinNumSpikes[n - 1], (spikes.array() != 0).count());
  }
}

}  // namespace audition
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
