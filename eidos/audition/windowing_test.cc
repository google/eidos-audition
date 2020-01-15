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

#include <cmath>
#include <vector>

#include "Eigen/Core"
#include "eidos/audition/stimulus_config.pb.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

using testing::Each;
using testing::Eq;

namespace eidos {
namespace audition {

constexpr int kNumChannels = 251;
constexpr int kNumSamples = 4000;
constexpr int kSampleRate = 16000;      // 16 kHz.
constexpr double kFrameShift = 1E-2;   // 10 ms.
constexpr double kWindowWidth = 25E-3;  // 25 ms.
constexpr int kNumExtraSamples = 3000;
constexpr int kNumExtraStep = 33;

TEST(WindowingTest, BasicCheckForFrames) {
  StimulusConfig config;
  config.set_sample_rate(kSampleRate);
  config.set_window_duration_sec(kWindowWidth);
  config.set_frame_shift_sec(kFrameShift);
  for (int i = 0; i < kNumExtraSamples; i += kNumExtraStep) {
    const int num_samples = kNumSamples + i;
    const Eigen::ArrayXXd input = Eigen::ArrayXXd::Random(
        kNumChannels, num_samples);
    const FrameInfo &info = GetFrameInfo(input, config);
    EXPECT_GT(info.num_frames, 0);
    EXPECT_GT(info.frame_size, 0);
    EXPECT_GT(info.frame_shift, 0);
    const std::vector<Eigen::ArrayXXd> frames = Window(input, config);
    EXPECT_EQ(frames.size(), info.num_frames);
    for (int j = 0; j < info.num_frames; ++j) {
      EXPECT_EQ(frames[j].rows(), kNumChannels);
      EXPECT_EQ(frames[j].cols(), info.frame_size);
    }
  }
}

TEST(WindowingTest, BasicTimeIntegration) {
  const Eigen::ArrayXXd input = Eigen::ArrayXXd::Random(
      kNumChannels, kNumSamples);
  StimulusConfig config;
  config.set_sample_rate(kSampleRate);
  config.set_window_duration_sec(kWindowWidth);
  config.set_frame_shift_sec(kFrameShift);
  const FrameInfo &info = GetFrameInfo(input, config);
  EXPECT_GT(info.num_frames, 0);
  const Eigen::ArrayXXd frames = WindowAndIntegrateTime(input, config);
  EXPECT_EQ(frames.cols(), info.num_frames);
  EXPECT_EQ(frames.rows(), input.rows());
}

TEST(WindowingTest, BasicWindowFunction) {
  constexpr int kFrameSize = 101;
  // No-op window function.
  Eigen::ArrayXXd window = ComputeWindowFunction(
      WINDOW_FUNCTION_NONE, kFrameSize, kNumChannels);
  EXPECT_EQ(window.rows(), kNumChannels);
  EXPECT_EQ(window.cols(), kFrameSize);
  const std::vector<double> window_vec(window.data(),
                                       window.data() + window.size());
  EXPECT_THAT(window_vec, Each(Eq(1.0)));
  // Hann window.
  window = ComputeWindowFunction(
      WINDOW_FUNCTION_HANN, kFrameSize, kNumChannels);
  EXPECT_EQ(window.rows(), kNumChannels);
  EXPECT_EQ(window.cols(), kFrameSize);
  EXPECT_EQ(1.0, window(0, std::ceil(kFrameSize / 2)));
}

}  // namespace audition
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
