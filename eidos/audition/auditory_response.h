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

// Simple handle encapsulating response from an auditory pipeline.

#ifndef EIDOS_AUDITION_AUDITORY_RESPONSE_H_
#define EIDOS_AUDITION_AUDITORY_RESPONSE_H_

#include <utility>
#include <string>
#include <vector>

#include "Eigen/Core"
#include "eidos/audition/auditory_model.h"
#include "eidos/audition/auditory_model_config.pb.h"

namespace eidos {
namespace audition {

class AuditoryResponse {
 public:
  AuditoryResponse() = default;
  ~AuditoryResponse() = default;

  // Writes the Eigen matrix to the supplied file in Numpy .npy format.
  void WriteNpy(const std::string &filename) const;

  // Writes a named collection of Eigen matrixes to the supplied file in Numpy
  // .npz format.
  void WriteNpz(const std::string &filename) const;

  const std::vector<std::pair<AuditoryOutputType, Eigen::ArrayXXd>> &
      outputs() const { return outputs_; }

  std::vector<std::pair<AuditoryOutputType, Eigen::ArrayXXd>> *
      mutable_outputs() { return &outputs_; }

  const ChannelProperties &channel_properties() const {
    return channel_properties_; }

  ChannelProperties *mutable_channel_properties() {
    return &channel_properties_; }

 private:
  // Collection of the outputs of all the stages in auditory chain.
  std::vector<std::pair<AuditoryOutputType, Eigen::ArrayXXd>> outputs_;

  // Encapsulation of various channel properties.
  ChannelProperties channel_properties_;
};

}  // namespace audition
}  // namespace eidos

#endif  // EIDOS_AUDITION_AUDITORY_RESPONSE_H_

// Local Variables:
// mode: c++
// End:
