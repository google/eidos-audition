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

// Simple abstraction representing an auditory pipeline consisting of
// auditory models.

#ifndef EIDOS_AUDITION_AUDITORY_PIPELINE_H_
#define EIDOS_AUDITION_AUDITORY_PIPELINE_H_

#include <memory>
#include <utility>
#include <vector>

#include "Eigen/Core"
#include "eidos/audition/auditory_model.h"
#include "eidos/audition/auditory_model_config.pb.h"
#include "eidos/audition/auditory_pipeline_config.pb.h"
#include "eidos/audition/auditory_response.h"
#include "eidos/audition/stimulus_config.pb.h"

namespace eidos {
namespace audition {

class AuditoryPipeline {
 public:
  AuditoryPipeline() = default;
  ~AuditoryPipeline() = default;

  // Initializes the pipeline from the configuration proto message.
  bool Init(const AuditoryPipelineConfig &config);

  // Initializes the pipeline given the list of models.
  bool Init(const std::vector<AuditoryModelType> &model_types);

  // Returns size of the pipeline corresponding to the number of models.
  int Size() const { return models_.size(); }

  // Feeds the stimulus to the pipeline. The stimulus described by the
  // configuration message <stimulus_config> is in <input>. Returns the output
  // in <response>. Depending on the settings, the response may contain a
  // single output from the last stage of auditory pipeline or outputs from all
  // the stages. The outputs always correspond to center frequencies arranged in
  // ascending order.
  void Process(const StimulusConfig &stimulus_config,
               const Eigen::Ref<const Eigen::ArrayXXd> &input,
               AuditoryResponse *response);

 protected:
  // Verifies the pipeline checking that models fit together.
  bool Verify() const;

 private:
  // Model pipeline.
  std::vector<std::unique_ptr<AuditoryModel>> models_;
};

}  // namespace audition
}  // namespace eidos

#endif  // EIDOS_AUDITION_AUDITORY_PIPELINE_H_

// Local Variables:
// mode: c++
// End:
