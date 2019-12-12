// Copyright 2000, 2010 Lowel P. O'Mard
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

// The hair cell synapse model from Sumner, et. al. (2002). The implementation
// is based on the reference implementation in DSAM known as <IHC_Meddis2000>
// model.
//
// References:
// -----------
// Sumner, C. J, Lopez-Poveda, E. A., O'Mard, L. P. and Meddis, R. (2002):
//   "A revised model of the inner-hair cell and auditory-nerve complex.",
//   The Journal of the Acoustical Society of America (JASA), vol.111, no.5,
//   pp. 2178--2188.

#ifndef THIRD_PARTY_AUDITION_MODELS_DSAM_SUMNER2002_HAIR_CELL_SYNAPSE_MODEL_H_
#define THIRD_PARTY_AUDITION_MODELS_DSAM_SUMNER2002_HAIR_CELL_SYNAPSE_MODEL_H_

#include "Eigen/Core"
#include "eidos/audition/auditory_model.h"
#include "eidos/audition/stimulus_config.pb.h"
#include "eidos/audition/sumner2002_hair_cell_synapse_config.pb.h"

namespace eidos {
namespace audition {

class Sumner2002HairCellSynapseModel : public AuditoryModel {
 public:
  explicit Sumner2002HairCellSynapseModel(
      const Sumner2002HairCellSynapseConfig &config);
  Sumner2002HairCellSynapseModel() = delete;
  ~Sumner2002HairCellSynapseModel() override = default;

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

 private:
  // Model configuration.
  Sumner2002HairCellSynapseConfig config_;
};

// Computes and returns the per-channel spike rate probabilities given the
// instantaneous amplitude of the signal <input> representing the basilar
// membrane displacement in each channel. The rows of the input stimulus
// correspond to channels along the cochlear partition, the columns to the
// response across time.
//
// The output has the same dimension as inputs. Each row in the output
// corresponds to the rate probability of neural spikes (spikes/sec)
// corresponding to each channel in the input.
void Sumner2002HairCellSynapse(
    const Sumner2002HairCellSynapseConfig &config,
    const Eigen::Ref<const Eigen::ArrayXXd> &input,
    int sample_rate,
    Eigen::ArrayXXd *output);

}  // namespace audition
}  // namespace eidos

#endif  // THIRD_PARTY_AUDITION_MODELS_DSAM_SUMNER2002_HAIR_CELL_SYNAPSE_MODEL_H_

// Local Variables:
// mode: c++
// End:
