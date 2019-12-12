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

// Interface for auditory models. Concrete filterbank implementations derive
// from this class.

#ifndef EIDOS_AUDITION_AUDITORY_MODEL_H_
#define EIDOS_AUDITION_AUDITORY_MODEL_H_

#include <string>
#include <vector>

#include "Eigen/Core"
#include "eidos/audition/auditory_model_config.pb.h"
#include "eidos/audition/stimulus_config.pb.h"

namespace Eigen {

// This definition seems to be missing from Eigen.
typedef Eigen::Array<double, 1, Eigen::Dynamic> RowArrayXd;

}  // namespace Eigen

namespace eidos {
namespace audition {

// Encapsulation of various channel properties.
struct ChannelProperties {
  // Returns number of channels.
  int NumChannels() const { return center_frequencies.size(); }

  // Center frequencies (CF) for each channel.
  std::vector<double> center_frequencies;
};

class AuditoryModel {
 public:
  AuditoryModel() = delete;
  virtual ~AuditoryModel() = default;

  // Initializes the model from the model-agnostic configuration.
  virtual void Init(const StimulusConfig &stimulus_config);

  // Returns the channel properties for the model.
  virtual ChannelProperties GetChannelProperties() const = 0;

  // When the input is audio, the audio is assumed to be monaural (single ear).
  // In this case the input has size (number of input channels or "ears", which
  // is 1) by (number of samples). When input is audio, it is assumed to be
  // normalized within [-1.0, 1.0] range and the gain already applied.
  //
  // Note that since the memory layout of an Eigen:ArrayXXd is
  // column-major, this means that all of the values for the first
  // time step are stored before all of the values for the second
  // time step, etc.
  //
  // The output is an ArrayXXd of size GetNumberOfFrequencyChannels() by
  // input.cols().
  virtual void ProcessSegment(const Eigen::Ref<const Eigen::ArrayXXd> &input,
                              Eigen::ArrayXXd *output) = 0;

  AuditoryModelType model_type() const { return model_type_; }
  AuditoryStageType stage_type() const { return stage_type_; }
  AuditoryOutputType output_type() const { return output_type_; }
  bool waveform_input() const { return waveform_input_; }
  const StimulusConfig &stimulus_config() const { return stimulus_config_; }
  StimulusConfig *mutable_stimulus_config() { return &stimulus_config_; }

  const AuditoryModel *previous_model() const { return previous_model_; }
  void set_previous_model(const AuditoryModel *model) {
    previous_model_ = model; }

  // Returns model type as string.
  std::string ModelTypeAsString() const;

 protected:
  AuditoryModel(AuditoryModelType model_type, AuditoryStageType stage_type,
                AuditoryOutputType output_type, bool waveform_input);

 private:
  // Model type.
  const AuditoryModelType model_type_;

  // Auditory stage that this model supports.
  const AuditoryStageType stage_type_;

  // Type of the output signal that model outputs. In theory we should be
  // supporting multiple outputs for a single model (especially for the
  // synapse), but keep it simple for now.
  const AuditoryOutputType output_type_;

  // If enabled, the input stimulus is the pressure waveform, i.e. this is the
  // head model. Otherwise, this is an intermediate model in the auditory chain.
  const bool waveform_input_;

  // Model-agnostic configuration.
  StimulusConfig stimulus_config_;

  // If models are in a cascade, keep the pointer to the previous model.
  // The pointer is owned elsewhere.
  const AuditoryModel *previous_model_;
};

// Returns auditory stage corresponding to the given auditory output type.
AuditoryStageType OutputType2Stage(AuditoryOutputType output_type);

}  // namespace audition
}  // namespace eidos

#endif  // EIDOS_AUDITION_AUDITORY_MODEL_H_

// Local Variables:
// mode: c++
// End:
