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

#include "eidos/audition/auditory_pipeline.h"

#include <algorithm>
#include <functional>
#include <map>
#include <set>
#include <string>

#include "absl/memory/memory.h"
#include "eidos/audition/carfac_config.pb.h"
#include "eidos/audition/carfac_model.h"
#include "eidos/audition/gammatone_filterbank_config.pb.h"
#include "eidos/audition/mock_basilar_membrane_model.h"
#include "eidos/audition/sumner2002_hair_cell_synapse_config.pb.h"
#include "eidos/audition/zhang2001_spike_generator_config.pb.h"
#include "eidos/stubs/integral_types.h"
#include "eidos/stubs/logging.h"
#include "third_party/audition/models/amt/meddis1986_hair_cell_synapse_model.h"
#include "third_party/audition/models/baumgarte/baumgarte_model.h"
#include "third_party/audition/models/bruce_carney/bruce2018_hair_cell_synapse_model.h"
#include "third_party/audition/models/bruce_carney/bruce2018_spike_generator.h"
#include "third_party/audition/models/bruce_carney/jackson_spike_generator.h"
#include "third_party/audition/models/bruce_carney/zhang2001_spike_generator.h"
#include "third_party/audition/models/bruce_carney/zilany2014_ihc_model.h"
#include "third_party/audition/models/dsam/sumner2002_hair_cell_synapse_model.h"
#include "third_party/audition/models/gammatone/gammatone_filterbank_model.h"

namespace eidos {
namespace audition {
namespace {

// Extracts model-specific configuration from the generic configuration message.
// If no configuration is present, the defaults will be passed to the model.
template <typename ModelConfigType>
ModelConfigType GetModelConfig(const AuditoryModelConfig &config) {
  ModelConfigType model_config;
  if (config.has_config() && config.config().Is<ModelConfigType>()) {
    const std::string model_name = AuditoryModelType_Name(config.model_type());
    GOOGLE_CHECK(config.config().UnpackTo(&model_config))
        << model_name << ": Failed to extract configuration";
    GOOGLE_LOG(INFO) << "Configuration [" << model_name << "]: "
                     << model_config.DebugString();
  }
  return model_config;
}

// Given the model configuration generates instance of the model.
std::unique_ptr<AuditoryModel> MakeModel(const AuditoryModelConfig &config) {
  const AuditoryModelType model_type = config.model_type();
  switch (model_type) {
    case MODEL_MOCK_BASILAR_MEMBRANE:
      return absl::make_unique<MockBasilarMembraneModel>();
    case MODEL_CARFAC:
      return absl::make_unique<CarfacModel>(
          GetModelConfig<CarfacConfig>(config));
    case MODEL_BAUMGARTE:
      return absl::make_unique<BaumgarteModel>();
    case MODEL_GAMMATONE_SLANEY:
      return absl::make_unique<GammatoneFilterbankModel>(
          GetModelConfig<GammatoneFilterbankConfig>(config));
    case MODEL_ZILANY_IHC_2014:
      return absl::make_unique<Zilany2014IhcModel>();
    case MODEL_MEDDIS_SYNAPSE_1986:
      return absl::make_unique<Meddis1986HairCellSynapseModel>();
    case MODEL_BRUCE_SYNAPSE_2018:
      return absl::make_unique<Bruce2018HairCellSynapseModel>();
    case MODEL_SUMNER_SYNAPSE_2002:
      return absl::make_unique<Sumner2002HairCellSynapseModel>(
          GetModelConfig<Sumner2002HairCellSynapseConfig>(config));
    case MODEL_BRUCE_SPIKES_2018:
      return absl::make_unique<Bruce2018SpikeGeneratorModel>();
    case MODEL_JACKSON_SPIKES:
      return absl::make_unique<JacksonSpikeGeneratorModel>();
    case MODEL_ZHANG_SPIKES_2001:
      return absl::make_unique<Zhang2001SpikeGeneratorModel>(
          GetModelConfig<Zhang2001SpikeGeneratorConfig>(config));
    default: {
      GOOGLE_LOG(FATAL) << "Unknown model type";
      return std::unique_ptr<AuditoryModel>(nullptr);
    }
  }
}

}  // namespace

bool AuditoryPipeline::Init(const AuditoryPipelineConfig &config) {
  if (config.models().empty()) {
    GOOGLE_LOG(ERROR) << "No models were supplied";
    return false;
  }
  models_.reserve(config.models().size());
  int model_counter = 0;
  std::set<AuditoryModelType> models_set;
  for (const auto &model_config : config.models()) {
    const AuditoryModelType model_type = model_config.model_type();
    if (model_type == MODEL_UNKNOWN) {
      GOOGLE_LOG(ERROR) << "Model " << model_counter << ": Unknown model type!";
      return false;
    }
    const std::string &model_name = AuditoryModelType_Name(model_type);
    if (models_set.find(model_type) != models_set.end()) {
      GOOGLE_LOG(ERROR) << "Model " << model_counter << ": Model " << model_name
                        << " is already in the pipeline";
      return false;
    }
    models_set.insert(model_type);
    std::unique_ptr<AuditoryModel> model = MakeModel(model_config);
    GOOGLE_LOG(INFO) << "Creating " << model_name << " ...";
    models_.emplace_back(std::move(model));
    model_counter++;
  }
  // Create the links between the models.
  for (uint32 i = 0; i < models_.size(); ++i) {
    if (i > 0) models_[i]->set_previous_model(models_[i - 1].get());
  }
  return Verify();
}

bool AuditoryPipeline::Init(const std::vector<AuditoryModelType> &model_types) {
  if (model_types.empty()) {
    GOOGLE_LOG(ERROR) << "No models were supplied";
    return false;
  }
  AuditoryPipelineConfig config;
  config.mutable_models()->Reserve(model_types.size());
  for (const auto model_type : model_types) {
    // Following will generate model with an empty config.
    AuditoryModelConfig *model_config = config.add_models();
    model_config->set_model_type(model_type);
  }
  return Init(config);
}

bool AuditoryPipeline::Verify() const {
  if (!models_[0]->waveform_input()) {
    GOOGLE_LOG(ERROR) << "The head model " << models_[0]->ModelTypeAsString()
                      << " does not take waveform as input";
    return false;
  }
  // A very basic check that the models can connect to each other provided by
  // the mapping between current auditory stage and a set of outputs that can
  // feed into it. This mechanism is not perfect, but will do for now.
  static const std::map<AuditoryStageType, std::set<AuditoryOutputType>>
      kStageValidLinks = {
    { STAGE_HAIR_CELLS, { OUTPUT_BASILAR_MEMBRANE_DISPLACEMENT } },
    { STAGE_AUDITORY_NERVE_SYNAPSE, { OUTPUT_BASILAR_MEMBRANE_DISPLACEMENT,
                                      OUTPUT_IHC_TRANSMEMBRANE_POTENTIAL,
                                      OUTPUT_SYNAPSE_FIRING_RATES} } };
  for (uint32 i = 1; i < models_.size(); ++i) {
    const AuditoryStageType stage_type = models_[i]->stage_type();
    const auto allowed_inputs_pos = kStageValidLinks.find(stage_type);
    if (allowed_inputs_pos == kStageValidLinks.end()) {
      GOOGLE_LOG(ERROR) << "Model " << models_[0]->ModelTypeAsString()
                        << ": No link constraints found for stage "
                        << AuditoryStageType_Name(stage_type);
      return false;
    }
    const AuditoryOutputType input_type = models_[i - 1]->output_type();
    const auto &allowed_inputs = allowed_inputs_pos->second;
    if (allowed_inputs.find(input_type) == allowed_inputs.end()) {
      GOOGLE_LOG(ERROR) << "Model " << models_[0]->ModelTypeAsString()
                        << ": Can't have "
                        << AuditoryOutputType_Name(input_type) << " as input";
      return false;
    }
  }
  GOOGLE_LOG(INFO) << "Loaded pipeline with " << models_.size() << " models.";
  return true;
}

void AuditoryPipeline::Process(const StimulusConfig &stimulus_config,
                               const Eigen::Ref<const Eigen::ArrayXXd> &input,
                               AuditoryResponse *response) {
  GOOGLE_CHECK(!models_.empty()) << "Pipeline is empty!";
  GOOGLE_CHECK(response != nullptr);
  auto &outputs = *response->mutable_outputs();
  outputs.clear();
  if (stimulus_config.store_intermediate_outputs()) {
    outputs.reserve(models_.size());
  } else {
    outputs.reserve(1);
  }
  Eigen::ArrayXXd model_input = input, model_output;
  const int num_models = models_.size();
  for (int i = 0; i < num_models; ++i) {
    AuditoryModel *model = models_[i].get();
    model->Init(stimulus_config);
    GOOGLE_LOG(INFO) << "Running " << model->ModelTypeAsString() << " ...";
    model->ProcessSegment(model_input, &model_output);
    if (stimulus_config.store_intermediate_outputs() || i == num_models - 1) {
      outputs.emplace_back(model->output_type(), model_output);
    }
    model_input.swap(model_output);
  }

  // Check if the channel frequencies sorted in descending order and always
  // present them in ascending order which may require rearranging the outputs.
  ChannelProperties channel_properties = models_.back()->GetChannelProperties();
  auto &center_frequencies = channel_properties.center_frequencies;
  if (std::is_sorted(center_frequencies.begin(), center_frequencies.end(),
                     std::greater<double>())) {  // Weak descending order.
    GOOGLE_LOG(INFO) << "Re-arranging CFs into ascending order ...";
    std::sort(center_frequencies.begin(), center_frequencies.end());
    const int num_channels = center_frequencies.size();
    int num_rearranged = 0;
    for (auto &model_output_info : outputs) {
      Eigen::ArrayXXd &model_outputs = model_output_info.second;
      if (model_outputs.rows() == num_channels) {
        // Re-arrange the output rows according to the order of CFs.
        model_outputs.colwise().reverseInPlace();
        num_rearranged++;
      }
    }
    GOOGLE_LOG(INFO) << num_rearranged << " outputs rearranged.";
  } else {
    GOOGLE_CHECK(std::is_sorted(center_frequencies.begin(),
                                center_frequencies.end()))
        << "The center frequencies are unordered!";
  }
  *response->mutable_channel_properties() = channel_properties;
}

}  // namespace audition
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
