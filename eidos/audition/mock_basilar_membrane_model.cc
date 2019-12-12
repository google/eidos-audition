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

#include "eidos/audition/mock_basilar_membrane_model.h"

#include "eidos/audition/psychoacoustic_scale_utils.h"

namespace eidos {
namespace audition {
namespace {

// Minimum characteristic frequency in Hz (human).
constexpr double kMinFrequency = 25.0;

// Maximum characteristic frequency in Hz (human).
constexpr double kMaxFrequency = 10E+3;

// Mock transmembrane voltages.
constexpr double kMinVoltage = -0.018051;  // ~18 mV
constexpr double kMaxVoltage = 0.038092;  // ~38 mV

// Creates mock properties for the model.
ChannelProperties CreateProperties(int num_channels) {
  ChannelProperties properties;
  properties.center_frequencies = erb::GetFrequenciesHz(
      num_channels, kMinFrequency, kMaxFrequency);
  return properties;
}

}  // namespace

constexpr int MockBasilarMembraneModel::kDefaultNumChannelsForWaveformInput;

MockBasilarMembraneModel::MockBasilarMembraneModel() :
    AuditoryModel(MODEL_MOCK_BASILAR_MEMBRANE, STAGE_BASILAR_MEMBRANE,
                  OUTPUT_BASILAR_MEMBRANE_DISPLACEMENT,
                  true /* waveform input */) {}

void MockBasilarMembraneModel::Init(const StimulusConfig &stimulus_config) {
  AuditoryModel::Init(stimulus_config);
}

ChannelProperties MockBasilarMembraneModel::GetChannelProperties() const {
  return channel_properties_;
}

void MockBasilarMembraneModel::ProcessSegment(
    const Eigen::Ref<const Eigen::ArrayXXd> &input, Eigen::ArrayXXd *output) {
  const int num_samples = input.cols();
  if (input.rows() == 1) {
    int num_channels = stimulus_config().num_channels();
    if (num_channels <= 0) {
      num_channels = kDefaultNumChannelsForWaveformInput;
    }
    channel_properties_ = CreateProperties(num_channels);
    // Simulate transmembrane voltages in [kMinVoltage, kMaxVoltage] V.
    *output = kMinVoltage + (Eigen::ArrayXXd::Random(
        num_channels, num_samples) * 0.5 + 0.5) *
        (kMaxVoltage - kMinVoltage);
  } else {
    channel_properties_ = CreateProperties(input.rows());
    *output = input;
  }
}

}  // namespace audition
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
