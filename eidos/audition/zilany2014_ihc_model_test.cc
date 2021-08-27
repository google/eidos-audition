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

// Simple test that the IHC model from Zilany, et. al. (2014) behaves.

#include <algorithm>

#include "Eigen/Core"
#include "eidos/audition/stimulus_config.pb.h"
#include "eidos/port/logging.h"
#include "gtest/gtest.h"
#include "third_party/audition/models/bruce_carney/zilany2014_ihc_model.h"

namespace eidos {
namespace audition {
namespace {

// Checks the model over the random monaural waveform.
void CheckRandomStimulus(int requested_num_channels,
                         int expected_num_channels,
                         int min_output_level) {
  StimulusConfig stimulus_config;
  stimulus_config.set_sample_rate(100000 /* 100 kHz */);
  stimulus_config.set_num_channels(requested_num_channels);

  const int kNumEars = 1;
  const int kNumSamples = 5000;
  const Eigen::ArrayXXd input = Eigen::ArrayXXd::Random(kNumEars, kNumSamples);

  Zilany2014IhcModel model;
  model.Init(stimulus_config);
  Eigen::ArrayXXd output;
  model.ProcessSegment(input, &output);
  EXPECT_EQ(output.rows(), expected_num_channels);
  EXPECT_EQ(output.cols(), kNumSamples);
  const ChannelProperties &channel_properties = model.GetChannelProperties();
  EXPECT_EQ(output.rows(), channel_properties.NumChannels());
  EXPECT_EQ(output.size(), channel_properties.NumChannels() * kNumSamples);
  const auto &center_frequencies = channel_properties.center_frequencies;
  EXPECT_TRUE(std::is_sorted(center_frequencies.begin(),
                             center_frequencies.end()));
  const float output_level = output.matrix().squaredNorm();
  GOOGLE_LOG(INFO) << "Output level: " << output_level;
  EXPECT_LT(min_output_level, output_level);
  EXPECT_LT(-1.5, output.minCoeff());  // Minimal potential.
  EXPECT_GT(2.0, output.maxCoeff());   // Maximal potential.
}

}  // namespace

TEST(Zilany2014IhcModelTest, RandomizedWaveformWithDefaultsTest) {
  // Mimimal permissable value for the amplification in the output.
  constexpr double kMinOutputLevel = 910.0;
  CheckRandomStimulus(0 /* requested channels */, 251 /* actual channels */,
                      kMinOutputLevel);
}

TEST(Zilany2014IhcModelTest, RandomizedWaveformTest) {
  constexpr int kNumChannels = 100;
  constexpr double kMinOutputLevel = 350.0;
  CheckRandomStimulus(kNumChannels /* requested channels */,
                      kNumChannels /* actual channels */,
                      kMinOutputLevel);
}

}  // namespace audition
}  // namespace eidos
