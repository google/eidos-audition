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

#include "eidos/audition/auditory_model.h"

namespace eidos {
namespace audition {

AuditoryStageType OutputType2Stage(AuditoryOutputType output_type) {
  if (output_type == OUTPUT_BASILAR_MEMBRANE_DISPLACEMENT) {
    return STAGE_BASILAR_MEMBRANE;
  } else if (output_type == OUTPUT_IHC_TRANSMEMBRANE_POTENTIAL) {
    return STAGE_HAIR_CELLS;
  } else {
    return STAGE_AUDITORY_NERVE_SYNAPSE;
  }
}

AuditoryModel::AuditoryModel(AuditoryModelType model_type,
                             AuditoryStageType stage_type,
                             AuditoryOutputType output_type,
                             bool waveform_input) :
    model_type_(model_type), stage_type_(stage_type), output_type_(output_type),
    waveform_input_(waveform_input), previous_model_(nullptr) {
}

void AuditoryModel::Init(const StimulusConfig &stimulus_config) {
  stimulus_config_ = stimulus_config;
}

std::string AuditoryModel::ModelTypeAsString() const {
  return AuditoryModelType_Name(model_type_);
}

}  // namespace audition
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
