// Copyright 2001 Xuedong Zhang
// Copyright 2019 Google LLC
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// version 2 as published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.

// Spike generation algorithm described in Zhang, X., Heinz, M. G., Bruce,
// I. C., & Carney, L. H. (2001): "A phenomenological model for the responses
// of auditory-nerve fibers: I. Nonlinear tuning with compression and
// suppression.", The Journal of the Acoustical Society of America,
// 109(2), 648-670.

#ifndef THIRD_PARTY_AUDITION_MODELS_BRUCE_CARNEY_ZHANG2001_SPIKE_GENERATOR_H_
#define THIRD_PARTY_AUDITION_MODELS_BRUCE_CARNEY_ZHANG2001_SPIKE_GENERATOR_H_

#include "Eigen/Core"
#include "eidos/audition/auditory_model.h"
#include "eidos/audition/stimulus_config.pb.h"
#include "eidos/audition/zhang2001_spike_generator_config.pb.h"

namespace eidos {
namespace audition {

class Zhang2001SpikeGeneratorModel : public AuditoryModel {
 public:
  explicit Zhang2001SpikeGeneratorModel(
      const Zhang2001SpikeGeneratorConfig &config);
  Zhang2001SpikeGeneratorModel() = delete;
  ~Zhang2001SpikeGeneratorModel() override = default;

  // Initializes the model from the model-agnostic configuration.
  void Init(const StimulusConfig &stimulus_config) override;

  // Returns the channel properties for the model.
  ChannelProperties GetChannelProperties() const override;

  // Processes an arbitrary number of samples corresponding to firing rate
  // probabilities, with input of size (number of input channels) by (number of
  // samples).
  //
  // The output is an ArrayXXd of size GetNumberOfFrequencyChannels() by
  // input.cols(). It stores the firing times as fire/no fire boolean mask.
  void ProcessSegment(const Eigen::Ref<const Eigen::ArrayXXd> &input,
                      Eigen::ArrayXXd *output) override;

 private:
  // Model configuration.
  Zhang2001SpikeGeneratorConfig config_;
};

// Processes an arbitrary number of samples corresponding to firing rate
// probabilities in <input>, with input of size (number of input channels) by
// (number of samples).
//
// The output is an ArrayXXd of size GetNumberOfFrequencyChannels() by
// input.cols(). It stores the firing times as fire/no fire boolean mask.
void Zhang2001SpikeGenerator(const Eigen::Ref<const Eigen::ArrayXXd> &input,
                             double sample_period,  // Seconds.
                             int num_repeats,  // Number of trials.
                             Eigen::ArrayXXd *output);

}  // namespace audition
}  // namespace eidos

#endif  // THIRD_PARTY_AUDITION_MODELS_BRUCE_CARNEY_ZHANG2001_SPIKE_GENERATOR_H_

// Local Variables:
// mode: c++
// End:
