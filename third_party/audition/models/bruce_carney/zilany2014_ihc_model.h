// Copyright 2013 Muhammad S. A. Zilany, Ian C. Bruce, Rasha A. Ibrahim,
//                Paul C. Nelson and Laurel H. Carney
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

// Auditory periphery model of Zilany, Bruce, et. al (2014).

/* This is Version 5.2 of the code for auditory periphery model of:

    Zilany, M.S.A., Bruce, I.C., Nelson, P.C., and Carney, L.H. (2009). "A
        phenomenological model of the synapse between the inner hair cell and
        auditory nerve: Long-term adaptation with power-law dynamics, "Journal
        of the Acoustical Society of America 126(5): 2390-2412.

   with the modifications and simulation options described in:

    Zilany, M.S.A., Bruce, I.C., Ibrahim, R.A., and Carney, L.H. (2013).
        "Improved parameters and expanded simulation options for a model of the
        auditory periphery," in Abstracts of the 36th ARO Midwinter Research
        Meeting.

   Humanization in this version includes:
   - Human middle-ear filter, based on the linear middle-ear circuit model of
     Pascal et al. (JASA 1998)
   - Human BM tuning, based on Shera et al. (PNAS 2002) or Glasberg & Moore
     (Hear. Res. 1990)
   - Human frequency-offset of control-path filter (i.e., cochlear amplifier
     mechanism), based on Greenwood (JASA 1990)

   The modifications to the BM tuning are described in:

        Ibrahim, R. A., and Bruce, I. C. (2010). "Effects of peripheral tuning
            on the auditory nerve's representation of speech envelope and
            temporal fine structure cues," in The Neurophysiological Bases of
            Auditory Perception, eds. E. A. Lopez-Poveda and A. R. Palmer and
            R. Meddis, Springer, NY, pp. 429-438.

   Please cite these papers if you publish any research
   results obtained with this code or any modified versions of this code.

   See the file readme.txt for details of compiling and running the model.

   %%% (c) M. S. Arefeen Zilany (msazilany@gmail.com), Ian C. Bruce
   (ibruce@ieee.org), Rasha A. Ibrahim, Paul C. Nelson, and Laurel H. Carney -
   November 2013 %%%
*/

#ifndef THIRD_PARTY_AUDITION_MODELS_BRUCE_CARNEY_ZILANY2014_IHC_MODEL_H_
#define THIRD_PARTY_AUDITION_MODELS_BRUCE_CARNEY_ZILANY2014_IHC_MODEL_H_

#include "Eigen/Core"
#include "eidos/audition/auditory_model.h"
#include "eidos/audition/stimulus_config.pb.h"

namespace eidos {
namespace audition {

class Zilany2014IhcModel : public AuditoryModel {
 public:
  Zilany2014IhcModel();
  ~Zilany2014IhcModel() override = default;

  // Initializes the model from the model-agnostic configuration.
  void Init(const StimulusConfig &stimulus_config) override;

  // Returns the channel properties for the model.
  ChannelProperties GetChannelProperties() const override;

  // Processes an arbitrary number of monaural audio samples, with input of size
  // (number of input channels or "ears", which is 1) by (number of samples).
  //
  // The output is an ArrayXXd of size GetNumberOfFrequencyChannels() by
  // input.cols(). It stores the estimates of IHC transmembrane potentials.
  void ProcessSegment(const Eigen::Ref<const Eigen::ArrayXXd> &input,
                      Eigen::ArrayXXd *output) override;

 private:
  // Encapsulation of various channel properties.
  ChannelProperties channel_properties_;
};

}  // namespace audition
}  // namespace eidos

#endif  // THIRD_PARTY_AUDITION_MODELS_BRUCE_CARNEY_ZILANY2014_IHC_MODEL_H_

// Local Variables:
// mode: c++
// End:
