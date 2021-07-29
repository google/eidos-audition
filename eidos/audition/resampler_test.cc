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

#include "eidos/audition/resampler.h"

#include <vector>

#include "Eigen/Core"
#include "eidos/audition/stimulus_config.pb.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

using testing::AllOf;
using testing::Each;
using testing::Gt;
using testing::Le;

namespace eidos {
namespace audition {

constexpr int kNumChannels = 50;
constexpr int kNumSamples = 200000;
constexpr int kSampleRate = 96000;      // 96 kHz.
constexpr int kUpsamplingFactor = 1;
constexpr int kDownsamplingFactor = 6;  // 16 kHz = kSampleRate / 6.

TEST(WindowingTest, BasicCheckDownsampleStep) {
  // Check no-op.
  StimulusConfig config;
  config.set_downsample_step(1);
  const Eigen::ArrayXXd input = Eigen::ArrayXXd::Random(
      kNumChannels, kNumSamples);
  Eigen::ArrayXXd output = DownstepResample(input, config);
  EXPECT_EQ(input.cols(), output.cols());
  EXPECT_EQ(input.rows(), output.rows());

  config.set_downsample_step(3);
  output = DownstepResample(input, config);
  EXPECT_EQ(input.cols() / 3, output.cols());
  EXPECT_EQ(input.rows(), output.rows());
}

TEST(WindowingTest, BasicCheckNormalDownsample) {
  // Prepare configuration.
  StimulusConfig config;
  config.set_sample_rate(kSampleRate);
  config.set_output_resample_up_factor(kUpsamplingFactor);
  config.set_output_resample_down_factor(kDownsamplingFactor);
  EXPECT_TRUE(ShouldResampleOutputs(config));

  // Call the downsampler.
  constexpr double kOriginalVal = 1.0;
  const Eigen::ArrayXXd input = Eigen::ArrayXXd::Constant(
      kNumChannels, kNumSamples, kOriginalVal);
  const Eigen::ArrayXXd output = Resample(input, config);
  EXPECT_EQ(kNumChannels, output.rows());
  const int num_output_samples = output.cols();
  EXPECT_EQ(num_output_samples, kNumSamples / kDownsamplingFactor + 1);

  // Check resampled values, these should not be 1.0 because of the filtering.
  const std::vector<double> resampled(output.data(),
                                      output.data() + output.size());
  EXPECT_THAT(resampled, Each(AllOf(Gt(kOriginalVal / 2),
                                    Le(kOriginalVal * 2))));
}

}  // namespace audition
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
