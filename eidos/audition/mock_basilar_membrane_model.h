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

// Mock basilar membrane model. Leaves the input stimulus unmodified if it
// doesn't look like a monoaural waveform. Otherwise, generates random outputs
// with a specified (or fixed) number of channels.

#ifndef EIDOS_AUDITION_MOCK_BASILAR_MEMBRANE_MODEL_H_
#define EIDOS_AUDITION_MOCK_BASILAR_MEMBRANE_MODEL_H_

#include "Eigen/Core"
#include "eidos/audition/auditory_model.h"
#include "eidos/audition/stimulus_config.pb.h"

namespace eidos {
namespace audition {

class MockBasilarMembraneModel : public AuditoryModel {
 public:
  MockBasilarMembraneModel();
  virtual ~MockBasilarMembraneModel() = default;

  // Initializes the model from the model-agnostic configuration.
  // Currently we use the internal defaults.
  void Init(const StimulusConfig &stimulus_config) override;

  // Returns the channel properties for the model.
  ChannelProperties GetChannelProperties() const override;

  // If the input is a waveform, generates random output with a fixed number of
  // channels. Otherwise, passes the input to the output unmodified.
  void ProcessSegment(const Eigen::Ref<const Eigen::ArrayXXd> &input,
                      Eigen::ArrayXXd *output) override;

  // Number of simulated channels.
  static constexpr int kDefaultNumChannelsForWaveformInput = 100;

 private:
  // Encapsulation of various channel properties.
  ChannelProperties channel_properties_;
};

}  // namespace audition
}  // namespace eidos

#endif  // EIDOS_AUDITION_MOCK_BASILAR_MEMBRANE_MODEL_H_

// Local Variables:
// mode: c++
// End:
