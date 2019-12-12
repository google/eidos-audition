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

// Simple test for synaptic model by Sumner, et. al. (2002).

#include "Eigen/Core"
#include "eidos/audition/gammatone_filterbank_config.pb.h"
#include "eidos/audition/sumner2002_hair_cell_synapse_config.pb.h"
#include "gtest/gtest.h"
#include "third_party/audition/models/dsam/sumner2002_hair_cell_synapse_model.h"
#include "third_party/audition/models/gammatone/gammatone_filterbank_model.h"

namespace eidos {
namespace audition {
namespace {

// Fills in and returns model configuration.
Sumner2002HairCellSynapseConfig GetConfig(bool output_spikes) {
  Sumner2002HairCellSynapseConfig config;
  config.set_ca_cond_mode(SUMNER2002_CA_COND_MODE_ORIGINAL);
  config.set_cleft_replenish_mode(SUMNER2002_CLEFT_REPLENISH_MODE_ORIGINAL);
  config.set_output_spikes(output_spikes);
  return config;
}

}  // namespace

TEST(Sumner2002HairCellSynapseModelTest, CheckRandomOriginalPaperProbs) {
  // This model is extremely sensitive to the input stimulus levels. Supplying
  // "wrong" stimulus intensity will result in numeric overflows and/or invalid
  // probability estimates.
  constexpr int kNumChannels = 10;
  constexpr int kNumSamples = 200000;
  constexpr int kSampleRate = 100000;  // Hz.
  constexpr double kInputScaleFactor = 100.0;
  const Eigen::ArrayXXd input =
      Eigen::ArrayXXd::Random(kNumChannels, kNumSamples) / kInputScaleFactor;
  Eigen::ArrayXXd output;
  const Sumner2002HairCellSynapseConfig config = GetConfig(
      /* output_spikes */false);
  Sumner2002HairCellSynapse(config, input, kSampleRate, &output);
  EXPECT_EQ(kNumChannels, output.rows());
  EXPECT_EQ(kNumSamples, output.cols());
  EXPECT_LE(0.0, output.minCoeff());
  EXPECT_GE(1.0, output.maxCoeff());
}

TEST(Sumner2002HairCellSynapseModelTest, CheckRandomOriginalPaperSpikes) {
  constexpr int kNumChannels = 10;
  constexpr int kNumSamples = 200000;
  constexpr int kSampleRate = 100000;  // Hz.
  constexpr double kInputScaleFactor = 5.0;
  const Eigen::ArrayXXd input =
      Eigen::ArrayXXd::Random(kNumChannels, kNumSamples) / kInputScaleFactor;
  Eigen::ArrayXXd output;
  const Sumner2002HairCellSynapseConfig config = GetConfig(
      /* output_spikes */true);
  Sumner2002HairCellSynapse(config, input, kSampleRate, &output);
  EXPECT_EQ(kNumChannels, output.rows());
  EXPECT_EQ(kNumSamples, output.cols());
  EXPECT_LE(0.0, output.minCoeff());
  EXPECT_GE(5.0, output.maxCoeff());
}

TEST(Sumner2002HairCellSynapseTest, GammatoneOverRandomStimulus) {
  constexpr int kNumChannels = 251;
  constexpr int kNumSamples = 10000;  // Long enough so we get a good average.
  constexpr double kInputScaleFactor = 5.0;
  const Eigen::ArrayXXd input =
      Eigen::ArrayXXd::Random(kNumChannels, kNumSamples) / kInputScaleFactor;
  constexpr int kSampleRate = 48000;  // Hz.

  // Estimate basilar membrane displacements.
  StimulusConfig stimulus_config;
  stimulus_config.set_sample_rate(kSampleRate /* Hz */);
  GammatoneFilterbankConfig gammatone_config;
  gammatone_config.set_filter_type(GAMMATONE_FILTERBANK_SLANEY);
  GammatoneFilterbankModel model(gammatone_config);
  model.Init(stimulus_config);
  Eigen::ArrayXXd bm_output;
  model.ProcessSegment(input, &bm_output);
  const ChannelProperties &channel_properties = model.GetChannelProperties();
  EXPECT_EQ(channel_properties.NumChannels(), bm_output.rows());
  EXPECT_EQ(kNumSamples, bm_output.cols());

  // Run Sumner's model.
  Eigen::ArrayXXd output;
  const Sumner2002HairCellSynapseConfig config = GetConfig(
      /* output_spikes */false);
  Sumner2002HairCellSynapse(config, bm_output, kSampleRate, &output);
  EXPECT_EQ(kNumChannels, output.rows());
  EXPECT_EQ(kNumSamples, output.cols());
  EXPECT_LE(0.0, output.minCoeff());
  EXPECT_GE(1.0, output.maxCoeff());
}

}  // namespace audition
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
