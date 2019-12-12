// Copyright 2007 Ian C. Bruce, Yousof Erfani, Muhammad S. A. Zilany
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

// Spike generator from the auditory nerve synapse model Carney, Bruce and
// Zilany labs. The model is based on the version described in:
//
//   Bruce, I.C., Erfani, Y., and Zilany, M.S.A. (2018). "A Phenomenological
//   model of the synapse between the inner hair cell and auditory nerve:
//   Implications of limited neurotransmitter release sites", Hearing research,
//   360, 40--54, (Special Issue on "Computational Models in Hearing").
//
// See the original documentation and copyright in bruce2018_model.README.

#ifndef THIRD_PARTY_AUDITION_MODELS_BRUCE_CARNEY_BRUCE2018_SPIKE_GENERATOR_H_
#define THIRD_PARTY_AUDITION_MODELS_BRUCE_CARNEY_BRUCE2018_SPIKE_GENERATOR_H_

#include <vector>

#include "Eigen/Core"
#include "eidos/audition/auditory_model.h"
#include "eidos/audition/stimulus_config.pb.h"

namespace eidos {
namespace audition {

class Bruce2018SpikeGeneratorModel : public AuditoryModel {
 public:
  Bruce2018SpikeGeneratorModel();
  ~Bruce2018SpikeGeneratorModel() override = default;

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
};

// Computes the spike times for a single AN fiber given the synapse output rate
// probabilities (in spikes/sec). The <sample_period> is reciprocal of the
// sampling rate. It defines the size of each bin.
//
// Note: Some parameters, such as spontaneous firing rate, are fiber-specific
// and depend on the type of the fiber population: low spontaneous rate (LSR),
// medium spontaneous rate (MSR) and high spontaneous rate (HSR).
std::vector<double> Bruce2018SpikeGenerator(
    const std::vector<double> &rates,  // Rate probabilities (spikes/sec).
    double sample_period,  // Sampling period of the rate function, in seconds.
    double spont_rate,     // Channel-specific spontaneous firing rate
                           // (spikes/sec).
    double absolute_refractory_period,   // In seconds.
    double relative_refractory_period);  // In seconds.

// Computes and returns per-channel spikes given the per-channel firing rates
// produced by the synapse model. The rows of <rates> correspond to channels
// along the cochlear partition, the columns to the firing rates across time.
// The sampling period of the rate function is given by <sample_period>.
//
// The output has the same dimension as inputs. Each row in the output
// corresponds to the binary indicator of whether a given channel fires or not
// at a given time. The single-channel API above is used to generate this
// information.
//
// TODO(agutkin): Spontaneous rate and refractory periods should be
//                channel-specific.
// TODO(agutkin): Sparse matrix for the outputs is more efficient.
void Bruce2018SpikeGenerator(
    const Eigen::Ref<const Eigen::ArrayXXd> &rates,
    double sample_period,  // Sampling period of the rate function, in seconds.
    double spont_rate,     // Channel-specific spontaneous firing rate
                           // (spikes/sec).
    double absolute_refractory_period,  // In seconds.
    double relative_refractory_period,  // In seconds.
    Eigen::ArrayXXd *output);

}  // namespace audition
}  // namespace eidos

#endif  // THIRD_PARTY_AUDITION_MODELS_BRUCE_CARNEY_BRUCE2018_SPIKE_GENERATOR_H_

// Local Variables:
// mode: c++
// End:
