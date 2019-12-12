// Copyright 1999,2000 Frank Baumgarte
// Copyright 2019 Google LLC
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

// Frank Baumgarte's model of hearing periphery described in detail in his
// PhD thesis:
//
//   F. Baumgarte: "Ein psychophysiologisches Gehoermodell zur Nachbildung
//   von Wahrnehmungsschwellen fuer die Audiocodierung",
//   (in German), PhD Dissertation, University of Hannover, 2000.
//
//   The peripheral ear model is based on the structure of
//   Eberhard Zwicker's "Analogmodell". That model consists
//   of analog electrical elements.
//
//   Wolfgang Peisl described the basic WDF algorithm in his
//   PhD thesis:
//   W. Peisl: "Beschreibung aktiver nichtlinearer Effekte der
//   peripheren Schallverarbeitung des Gehoers durch ein Rechnermodell", (in
//   German), Dissertation, Technical University of Munich.
//
// This code was originally available from Frank's page at University of
// Hannover but the link has been long dead. Please see the original
// documentation in "README.txt".

#ifndef THIRD_PARTY_AUDITION_MODELS_BAUMGARTE_BAUMGARTE_MODEL_H_
#define THIRD_PARTY_AUDITION_MODELS_BAUMGARTE_BAUMGARTE_MODEL_H_

#include <memory>

#include "Eigen/Core"
#include "eidos/audition/auditory_model.h"
#include "eidos/audition/stimulus_config.pb.h"

namespace eidos {
namespace audition {

struct BaumgarteHandle;

class BaumgarteModel : public AuditoryModel {
 public:
  BaumgarteModel();
  ~BaumgarteModel() override;

  // Initializes the model from the model-agnostic configuration.
  // Currently we use the internal defaults.
  void Init(const StimulusConfig &stimulus_config) override;

  // Returns the channel properties for the model.
  ChannelProperties GetChannelProperties() const override;

  // Processes an arbitrary number of monaural audio samples, with input of size
  // (number of input channels or "ears", which is 1) by (number of samples).
  // Note that since the memory layout of an Eigen:ArrayXXd is
  // column-major, this means that all of the values for the first
  // time step are stored before all of the values for the second
  // time step, etc.
  //
  // The output is an ArrayXXd of size GetNumberOfFrequencyChannels() by
  // input.cols(). It stores the estimates of IHC trans-membrane potentials.
  void ProcessSegment(const Eigen::Ref<const Eigen::ArrayXXd> &input,
                      Eigen::ArrayXXd *output) override;

 private:
  // Encapsulation of various channel properties.
  ChannelProperties channel_properties_;

  // Filterbank implementation.
  std::unique_ptr<BaumgarteHandle> impl_;
};

}  // namespace audition
}  // namespace eidos

#endif  // THIRD_PARTY_AUDITION_MODELS_BAUMGARTE_BAUMGARTE_MODEL_H_

// Local Variables:
// mode: c++
// End:
