// Copyright 1998 Malcolm Slaney (Interval Research Corporation)
// Copyright 2016 Piotr Majdak
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

// Inner hair cell synapse model by Ray Meddis, et. al.
//
// References:
// -----------
//   - Ray Meddis, Michael J. Hewitt, and Trevor M. Shackleton: "Implementation
//   details of a computation model of the inner hair‐cell auditory‐nerve
//   synapse", The Journal of the Acoustical Society of America 87, 1813 (1990).
//   - Ray Meddis (1986): "Simulation of mechanical to neural transduction in
//   the auditory receptor", Journal of the Acoustical Society of America 79(3),
//   702--711.
//
// Implementation is based on Malcolm Slaney's Meddis model in AuditoryToolbox
// and reimplementation in AMT. A slightly more elaborate implementation can be
// found in <MoIHC86Meddis> module of DSAM.

#ifndef THIRD_PARTY_AUDITION_MODELS_AMT_MEDDIS1986_HAIR_CELL_SYNAPSE_MODEL_H_
#define THIRD_PARTY_AUDITION_MODELS_AMT_MEDDIS1986_HAIR_CELL_SYNAPSE_MODEL_H_

#include "Eigen/Core"
#include "eidos/audition/auditory_model.h"
#include "eidos/audition/stimulus_config.pb.h"

namespace eidos {
namespace audition {

class Meddis1986HairCellSynapseModel : public AuditoryModel {
 public:
  Meddis1986HairCellSynapseModel();
  ~Meddis1986HairCellSynapseModel() override = default;

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

// Computes and returns the per-channel spike rate probabilities given the
// instantaneous amplitude of the signal <input> representing the basilar
// membrane displacement in each channel. The rows of the input stimulus
// correspond to channels along the cochlear partition, the columns to the
// response across time.
//
// Originally the Meddis model was intended to be used in tandem with a
// Gammatone filterbank (proposed by Roy Patterson, et. al.), but any model
// giving the estimates of basilar membrane displacements in response to
// the input stimulus can be used as well.
//
// The output has the same dimension as inputs. Each row in the output
// corresponds to the rate probability of neural spikes (spikes/sec)
// corresponding to each channel in the input.
//
// Generic (not necessarily related to the *instantaneous* firing rate)
// note on "rate probabilities":
//   "The rate is simply the expected number of events per some (time) unit
//    That could easily be larger than one, there is indeed no upper limit
//    on a rate, just make the time interval larger, then the rate becomes
//    larger."
void Meddis1986HairCellSynapse(const Eigen::Ref<const Eigen::ArrayXXd> &input,
                               int sample_rate, bool subtract_spontaneous_rate,
                               Eigen::ArrayXXd *output);

}  // namespace audition
}  // namespace eidos

#endif  // THIRD_PARTY_AUDITION_MODELS_AMT_MEDDIS1986_HAIR_CELL_SYNAPSE_MODEL_H_

// Local Variables:
// mode: c++
// End:
