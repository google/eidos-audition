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

// Inner hair cell synaptic model from Carney, Bruce and Zilany labs. The model
// is based on the version described in:
//
//   Bruce, I.C., Erfani, Y., and Zilany, M.S.A. (2018). "A phenomenological
//   model of the synapse between the inner hair cell and auditory nerve:
//   Implications of limited neurotransmitter release sites", Hearing research,
//   360, 40--54, (Special Issue on "Computational Models in Hearing").
//
// This library implements the transduction stage between the IHC transmembrane
// potentials and AN firing rates.
//
// See the original documentation and copyright in bruce2018_model.README.

#ifndef THIRD_PARTY_AUDITION_MODELS_BRUCE_CARNEY_BRUCE2018_HAIR_CELL_SYNAPSE_MODEL_H_
#define THIRD_PARTY_AUDITION_MODELS_BRUCE_CARNEY_BRUCE2018_HAIR_CELL_SYNAPSE_MODEL_H_

#include "Eigen/Core"
#include "eidos/audition/auditory_model.h"
#include "eidos/audition/stimulus_config.pb.h"

namespace eidos {
namespace audition {

class Bruce2018HairCellSynapseModel : public AuditoryModel {
 public:
  Bruce2018HairCellSynapseModel();
  ~Bruce2018HairCellSynapseModel() override = default;

  // Initializes the model from the model-agnostic configuration.
  void Init(const StimulusConfig &stimulus_config) override;

  // Returns the channel properties for the model.
  ChannelProperties GetChannelProperties() const override;

  // Processes an arbitrary number of samples corresponding to IHC transmembrane
  // potentials, with input of size (number of input channels) by (number of
  // samples).
  //
  // The output is an ArrayXXd of size GetNumberOfFrequencyChannels() by
  // input.cols(). It stores the firing rate probabilities of individual fibers.
  void ProcessSegment(const Eigen::Ref<const Eigen::ArrayXXd> &input,
                      Eigen::ArrayXXd *output) override;
};

// Computes estimates for the firing rate probabilities at each time step given
// the transmembrane potentials in <ihc_out>. The sample period of the stimulus
// is given by <sample_period>. The characteristic frequency of the fiber in Hz
// is specified by <center_freq>, <spont_rate> is the spontaneous firing rate
// (Hz). Whether to use exact (slow) or approximate (fast) implementation of the
// Power-Law function is controlled by <exact_power_law>.
void Bruce2018HairCellSynapse(
    const Eigen::RowArrayXd &ihc_out, double sample_period,
    double center_freq, double spont_rate,
    bool exact_power_law, Eigen::RowArrayXd *rates);

}  // namespace audition
}  // namespace eidos

#endif  // THIRD_PARTY_AUDITION_MODELS_BRUCE_CARNEY_BRUCE2018_HAIR_CELL_SYNAPSE_MODEL_H_

// Local Variables:
// mode: c++
// End:
