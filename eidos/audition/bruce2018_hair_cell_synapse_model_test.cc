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

// Simple test that the synaptic model from Bruce, et. al. (2018) behaves.

#include "Eigen/Core"
#include "eidos/audition/mock_basilar_membrane_model.h"
#include "eidos/audition/stimulus_config.pb.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "third_party/audition/models/baumgarte/baumgarte_model.h"
#include "third_party/audition/models/bruce_carney/bruce2018_hair_cell_synapse_model.h"

namespace eidos {
namespace audition {

// Transmembrane voltages.
constexpr double kMinVoltage = -0.018051;  // ~18 mV
constexpr double kMaxVoltage = 0.038092;  // ~38 mV

// Sanity checks for a single fiber given a random voltages stimulus using
// either exact or approximate Power-Law implementation.
void CheckSingleFiberOnRandomVoltages(bool exact_power_law,
                                      Eigen::RowArrayXd *rates) {
  constexpr int kNumSamples = 200000;  // Long enough so we get a good average.
  constexpr int kSampleRate = 100000;  // Hz.
  constexpr double kSpontRate = 60.0;  // Spontaneous firing rate (spikes/sec).
  constexpr double kCenterFreq = 1000.0;  // CF in Hz.

  // Simulate transmembrane voltages in the range [kMinVoltage, kMaxVoltage] V.
  const Eigen::RowArrayXd ihc_voltages =
      kMinVoltage + (Eigen::RowArrayXd::Random(kNumSamples) * 0.5 + 0.5) *
      (kMaxVoltage - kMinVoltage);
  Bruce2018HairCellSynapse(ihc_voltages, 1.0 / kSampleRate, kCenterFreq,
                           kSpontRate, exact_power_law, rates);
  EXPECT_EQ(kNumSamples, rates->size());
}

TEST(Bruce2018HairCellSynapseModelTest, RandomSingleFiberExactPowerLaw) {
  Eigen::RowArrayXd rates;
  CheckSingleFiberOnRandomVoltages(true /* exact Power-Law */, &rates);
  EXPECT_LE(0.0, rates.minCoeff());
  EXPECT_GT(5300.0, rates.maxCoeff());  // Hz.
}

TEST(Bruce2018HairCellSynapseModelTest, RandomSingleFiberApproxPowerLaw) {
  Eigen::RowArrayXd rates;
  CheckSingleFiberOnRandomVoltages(false /* exact Power-Law */, &rates);
  EXPECT_LE(0.0, rates.minCoeff());
  EXPECT_GT(5300.0, rates.maxCoeff());  // Hz.
}

TEST(Bruce2018HairCellSynapseModelTest, RandomAllFibers) {
  // Simulate transmembrane voltages and pass them unmodified using the mock
  // basilar membrane model.
  constexpr int kNumChannels = 40;
  constexpr int kNumSamples = 100000;
  const Eigen::ArrayXXd stimulus =
      kMinVoltage + (Eigen::ArrayXXd::Random(kNumChannels, kNumSamples) * 0.5 +
                     0.5) * (kMaxVoltage - kMinVoltage);

  constexpr int kSampleRate = 100000;  // Hz.
  StimulusConfig stimulus_config;
  stimulus_config.set_sample_rate(kSampleRate);
  MockBasilarMembraneModel bm_model;
  bm_model.Init(stimulus_config);
  Eigen::ArrayXXd ihc_voltages;  // IHC voltages mocked above.
  bm_model.ProcessSegment(stimulus, &ihc_voltages);

  // Run the model.
  Bruce2018HairCellSynapseModel model;
  model.set_previous_model(&bm_model);
  model.Init(stimulus_config);
  Eigen::ArrayXXd rates;  // firing rates.
  model.ProcessSegment(ihc_voltages, &rates);

  const ChannelProperties &channel_properties = model.GetChannelProperties();
  EXPECT_EQ(kNumChannels, channel_properties.NumChannels());
  EXPECT_EQ(kNumChannels, rates.rows());
  EXPECT_EQ(kNumSamples, rates.cols());
  EXPECT_LE(0.0, rates.minCoeff());
  EXPECT_GT(70000.0, rates.maxCoeff());  // Hz.
}

// Simulate the transmembrane voltages using Baumgarte model.
TEST(Bruce2018HairCellSynapseModelTest, BaumgarteOverRandomInputs) {
  // Random input pressure wave.
  const int kNumEars = 1;
  const int kNumSamples = 20000;  // Long enough so we get a good average.
  const Eigen::ArrayXXd input = Eigen::ArrayXXd::Random(kNumEars, kNumSamples);

  // Estimate IHC voltages.
  const int kSampleRate = 100000;  // Hz
  const int kNumChannels = 251;
  StimulusConfig stimulus_config;
  stimulus_config.set_sample_rate(kSampleRate);
  stimulus_config.set_downsample_step(1);  // No downsampling.
  BaumgarteModel bm_model;
  bm_model.Init(stimulus_config);
  Eigen::ArrayXXd voltages;
  bm_model.ProcessSegment(input, &voltages);
  EXPECT_EQ(kNumChannels, voltages.rows());
  EXPECT_EQ(kNumSamples, voltages.cols());
  EXPECT_LT(-1.0, voltages.minCoeff());
  EXPECT_GT(2.0, voltages.maxCoeff());

  // Run the synapse model.
  Bruce2018HairCellSynapseModel synapse_model;
  synapse_model.set_previous_model(&bm_model);
  synapse_model.Init(stimulus_config);
  Eigen::ArrayXXd rates;  // firing rates.
  synapse_model.ProcessSegment(voltages, &rates);

  const ChannelProperties &channel_properties =
      synapse_model.GetChannelProperties();
  EXPECT_EQ(kNumChannels, channel_properties.NumChannels());
  EXPECT_EQ(kNumChannels, rates.rows());
  EXPECT_EQ(kNumSamples, rates.cols());
  EXPECT_LE(0.0, rates.minCoeff());
  EXPECT_GT(7000.0, rates.maxCoeff());  // Hz.
}

}  // namespace audition
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
