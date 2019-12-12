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

// Simple test for Ray Meddis' model of inner hair cell synapse.

#include "Eigen/Core"
#include "eidos/audition/gammatone_filterbank_config.pb.h"
#include "eidos/audition/stimulus_config.pb.h"
#include "gtest/gtest.h"
#include "third_party/audition/models/amt/meddis1986_hair_cell_synapse_model.h"
#include "third_party/audition/models/baumgarte/baumgarte_model.h"
#include "third_party/audition/models/gammatone/gammatone_filterbank_model.h"

namespace eidos {
namespace audition {

// Computes IHC synapse response and performs various checks depending on
// whether the steady-state response is removed from the response or not.
void RandomStimulus(bool subtract_spontaneous_rate) {
  constexpr int kNumChannels = 251;
  constexpr int kNumSamples = 10000;  // Long enough so we get a good average.
  const Eigen::ArrayXXd input =
      Eigen::ArrayXXd::Random(kNumChannels, kNumSamples);
  constexpr int kSampleRate = 16000;  // Hz.
  Eigen::ArrayXXd output;
  Meddis1986HairCellSynapse(input, kSampleRate, subtract_spontaneous_rate,
                            &output);
  EXPECT_EQ(kNumChannels, output.rows());
  EXPECT_EQ(kNumSamples, output.cols());
  // Check the bounds on steady-state spontaneous firing rate and the maximum
  // firing rate. Spontaneous rate is so small that it doesn't make a huge
  // effect on the firing rate limits.
  EXPECT_LT(10.0, output.minCoeff());
  EXPECT_GT(80.0, output.maxCoeff());
}

TEST(Meddis1986HairCellSynapseTest, RandomizedStimulusWithSpontaneousRateTest) {
  RandomStimulus(false /* subtract_spontaneous_rate */);
}

TEST(Meddis1986HairCellSynapseTest,
     RandomizedStimulusWithNoSpontaneousRateTest) {
  RandomStimulus(true /* subtract_spontaneous_rate */);
}

TEST(Meddis1986HairCellSynapseTest, GammatoneOverRandomStimulus) {
  constexpr int kNumChannels = 251;
  constexpr int kNumSamples = 10000;  // Long enough so we get a good average.
  const Eigen::ArrayXXd input =
      Eigen::ArrayXXd::Random(kNumChannels, kNumSamples);
  constexpr int kSampleRate = 16000;  // Hz.

  // Estimate basilar membrane displacements.
  StimulusConfig stimulus_config;
  stimulus_config.set_sample_rate(kSampleRate /* Hz */);
  GammatoneFilterbankConfig config;
  config.set_filter_type(GAMMATONE_FILTERBANK_SLANEY);
  GammatoneFilterbankModel model(config);
  model.Init(stimulus_config);
  Eigen::ArrayXXd bm_output;
  model.ProcessSegment(input, &bm_output);
  const ChannelProperties &channel_properties = model.GetChannelProperties();
  EXPECT_EQ(channel_properties.NumChannels(), bm_output.rows());
  EXPECT_EQ(kNumSamples, bm_output.cols());

  // Run Meddis' model.
  Eigen::ArrayXXd output;
  Meddis1986HairCellSynapse(bm_output, kSampleRate,
                        false /* subtract_spontaneous_rate */, &output);
  EXPECT_EQ(kNumChannels, output.rows());
  EXPECT_EQ(kNumSamples, output.cols());
  EXPECT_LT(60.0, output.minCoeff());
  EXPECT_GT(70.0, output.maxCoeff());
}

TEST(Meddis1986HairCellSynapseTest, BaumgarteOverRandomStimulus) {
  // Random input pressure wave.
  const int kNumEars = 1;
  const int kNumSamples = 20000;  // Long enough so we get a good average.
  const Eigen::ArrayXXd input = Eigen::ArrayXXd::Random(kNumEars, kNumSamples);

  // Estimate IHC voltages. For the Meddis model, we ideally need the
  // BM displacements instead.
  const int kSampleRate = 100000;  // Hz
  const int kNumChannels = 251;
  StimulusConfig stimulus_config;
  stimulus_config.set_sample_rate(kSampleRate);
  stimulus_config.set_downsample_step(1);  // No downsampling.
  BaumgarteModel model;
  model.Init(stimulus_config);
  Eigen::ArrayXXd voltages;
  model.ProcessSegment(input, &voltages);
  EXPECT_EQ(kNumChannels, voltages.rows());
  EXPECT_EQ(kNumSamples, voltages.cols());
  EXPECT_LT(-1.0, voltages.minCoeff());
  EXPECT_GT(2.0, voltages.maxCoeff());

  // Run the synapse model. Perform fake voltage -> BM displacement scaling.
  Meddis1986HairCellSynapseModel synapse_model;
  synapse_model.Init(stimulus_config);
  Eigen::ArrayXXd rates;  // firing rates.
  constexpr double kVoltageToDisplacement = 2000.0;
  const Eigen::ArrayXXd bm_displacements = voltages * kVoltageToDisplacement;
  synapse_model.ProcessSegment(bm_displacements, &rates);

  const ChannelProperties &channel_properties = model.GetChannelProperties();
  EXPECT_EQ(kNumChannels, channel_properties.NumChannels());
  EXPECT_EQ(kNumChannels, rates.rows());
  EXPECT_EQ(kNumSamples, rates.cols());
  EXPECT_LT(1.0, rates.minCoeff());  // Hz.
  EXPECT_GT(210.0, rates.maxCoeff());  // Hz.
}

}  // namespace audition
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
