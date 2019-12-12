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

// Wrapper around CARFAC library.
//
// References:
// -----------
// [1] Lyon, R. (2017): "The CARFAC Digital Cochlear Model.", In "Human and
// Machine Hearing: Extracting Meaning from Sound", Chapter 15, (pp. 293-298),
// Cambridge: Cambridge University Press. doi:10.1017/9781139051699.020.
// [2] Lyon, R. (2011): "Using a Cascade of Asymmetric Resonators with
// Fast-Acting Compression as a Cochlear Model for Machine-Hearing
// Applications", Autumn Meeting of the Acoustical Society of Japan (2011),
// pp. 509-512.
// [3] Lyon, R. (2017): "Human and Machine Hearing", Cambridge University Press.

#ifndef EIDOS_AUDITION_CARFAC_MODEL_H_
#define EIDOS_AUDITION_CARFAC_MODEL_H_

#include <memory>

#include "Eigen/Core"
#include "eidos/audition/auditory_model.h"
#include "eidos/audition/carfac_config.pb.h"
#include "eidos/audition/stimulus_config.pb.h"

namespace eidos {
namespace audition {

struct CarfacHandle;

class CarfacModel : public AuditoryModel {
 public:
  explicit CarfacModel(const CarfacConfig &config);
  CarfacModel() = delete;
  ~CarfacModel() override;

  // Initializes the model from the model-agnostic configuration.
  // Currently we use the internal CARFAC defaults.
  void Init(const StimulusConfig &stimulus_config) override;

  // Returns the channel properties for the model.
  ChannelProperties GetChannelProperties() const override;

  // Processes an arbitrary number of samples, with input of size
  // (number of input channels or "ears") by (number of samples).
  // Note that since the memory layout of an Eigen:ArrayXXd is
  // column-major, this means that all of the values for the first
  // time step are stored before all of the values for the second
  // time step, etc.
  //
  // The input is assumed to be normalized within [-1.0, 1.0] range and the
  // gain already applied. The output is an ArrayXXd of size
  // GetNumberOfFrequencyChannels() by input.cols(). Depending on the type of
  // the output given in the configuration (please see <carfac_config.proto>),
  // the type of the returned output will be *one* of:
  //   - Neural activity patterns (NAPs),
  //   - Basilar membrane displacements,
  //   - Outer hair cell potentials,
  //   - Adaptive gain control levels.
  //
  // This is a limitation of the current implementation. Ideally multiple output
  // types should be returned, if configured.
  void ProcessSegment(const Eigen::Ref<const Eigen::ArrayXXd> &input,
                      Eigen::ArrayXXd *output) override;

  // Resets internal state of the filterbank. Call to "forget" the
  // effect of any previous input on the filterbank.
  void Reset();

 private:
  // Model configuration.
  CarfacConfig config_;

  // Encapsulation of various channel properties.
  ChannelProperties channel_properties_;

  // Actual model implementation.
  std::unique_ptr<CarfacHandle> impl_;
};

}  // namespace audition
}  // namespace eidos

#endif  // EIDOS_AUDITION_CARFAC_MODEL_H_

// Local Variables:
// mode: c++
// End:
