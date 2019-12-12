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

// Simple test for the mock basilar membrane model.

#include "Eigen/Core"
#include "eidos/audition/mock_basilar_membrane_model.h"
#include "eidos/audition/stimulus_config.pb.h"
#include "gtest/gtest.h"

namespace eidos {
namespace audition {
namespace {

TEST(MockBasilarMembraneModelTest, CheckTwoDimensionalStimulus) {
  StimulusConfig stimulus_config;
  stimulus_config.set_sample_rate(100000 /* Hz */);
  stimulus_config.set_audio_scaling_gain(1.0);
  constexpr int kNumChannels = 100;
  constexpr int kNumSamples = 10000;
  const Eigen::ArrayXXd input = Eigen::ArrayXXd::Random(
      kNumChannels, kNumSamples);
  MockBasilarMembraneModel model;
  model.Init(stimulus_config);
  Eigen::ArrayXXd output;
  model.ProcessSegment(input, &output);
  EXPECT_EQ(kNumChannels, output.rows());
  EXPECT_EQ(kNumSamples, output.cols());
}

// Checks the case when the input is assumed to be monaural input.
void CheckWaveformWithSimulatedChannels(int input_num_channels,
                                        int expected_num_channels) {
  StimulusConfig stimulus_config;
  stimulus_config.set_sample_rate(100000 /* Hz */);
  stimulus_config.set_audio_scaling_gain(1.0);
  stimulus_config.set_num_channels(input_num_channels);
  constexpr int kNumSamples = 10000;
  const Eigen::ArrayXXd input = Eigen::ArrayXXd::Random(
      1, kNumSamples);
  MockBasilarMembraneModel model;
  model.Init(stimulus_config);
  Eigen::ArrayXXd output;
  model.ProcessSegment(input, &output);
  EXPECT_EQ(expected_num_channels, output.rows());
  EXPECT_EQ(kNumSamples, output.cols());
}

TEST(MockBasilarMembraneModelTest, CheckMockWaveformWithDefaultNumChannels) {
  CheckWaveformWithSimulatedChannels(
      0, MockBasilarMembraneModel::kDefaultNumChannelsForWaveformInput);
}

TEST(MockBasilarMembraneModelTest, CheckMockWaveformWithConfiguredNumChannels) {
  constexpr int kNumChannels = 50;
  CheckWaveformWithSimulatedChannels(kNumChannels, kNumChannels);
}

}  // namespace
}  // namespace audition
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
