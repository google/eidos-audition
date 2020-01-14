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

#include "eidos/audition/windowing.h"

#include <vector>

#include "Eigen/Core"
#include "gtest/gtest.h"

namespace eidos {
namespace audition {

constexpr int kNumChannels = 251;
constexpr int kNumSamples = 4000;
constexpr int kSampleRate = 16000;      // 16 kHz.
constexpr double kFrameShift = 1E-2;   // 10 ms.
constexpr double kWindowWidth = 25E-3;  // 25 ms.
constexpr int kNumExtraSamples = 3000;
constexpr int kNumExtraStep = 33;

TEST(WindowingTest, BasicCheck) {
  for (int i = 0; i < kNumExtraSamples; i += kNumExtraStep) {
    const int num_samples = kNumSamples + i;
    const Eigen::ArrayXXd input = Eigen::ArrayXXd::Random(
        kNumChannels, num_samples);
    const FrameInfo &info = GetFrameInfo(input, kSampleRate, kWindowWidth,
                                         kFrameShift);
    EXPECT_GT(info.num_frames, 0);
    EXPECT_GT(info.frame_size, 0);
    EXPECT_GT(info.frame_shift, 0);
    const std::vector<Eigen::ArrayXXd> frames = Window(
        input, kSampleRate, kWindowWidth, kFrameShift);
    EXPECT_EQ(frames.size(), info.num_frames);
    for (int j = 0; j < info.num_frames; ++j) {
      EXPECT_EQ(frames[j].rows(), kNumChannels);
      EXPECT_EQ(frames[j].cols(), info.frame_size);
    }
  }
}

}  // namespace audition
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
