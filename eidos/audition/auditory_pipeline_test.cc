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

// Check that the auditory pipeline interface behaves.

#include <algorithm>
#include <functional>

#include "Eigen/Core"
#include "eidos/audition/auditory_model.h"
#include "eidos/audition/auditory_model_config.pb.h"
#include "eidos/audition/auditory_pipeline.h"
#include "eidos/audition/auditory_response.h"
#include "eidos/audition/stimulus_config.pb.h"
#include "gtest/gtest.h"

namespace eidos {
namespace audition {

TEST(AuditoryPipelineTest, CheckEmpty) {
  AuditoryPipeline pipeline;
  std::vector<AuditoryModelType> empty_models;
  EXPECT_FALSE(pipeline.Init(empty_models));
  AuditoryPipelineConfig empty_config;
  EXPECT_FALSE(pipeline.Init(empty_config));
  empty_config.add_models();  // Add empty model.
  EXPECT_FALSE(pipeline.Init(empty_config));
}

// Duplicates are not allowed.
TEST(AuditoryPipelineTest, CheckDuplicate) {
  AuditoryPipeline pipeline;
  const std::vector<AuditoryModelType> models = {
    MODEL_GAMMATONE_SLANEY, MODEL_GAMMATONE_SLANEY};
  EXPECT_FALSE(pipeline.Init(models));
}

// Checks that all head models take audio as an input.
TEST(AuditoryPipelineTest, CheckNoWaveformInputs) {
  const std::vector<AuditoryModelType> non_waveform_models = {
    MODEL_MEDDIS_SYNAPSE_1986, MODEL_JACKSON_SPIKES};
  for (const auto model_type : non_waveform_models) {
    std::vector<AuditoryModelType> models = { model_type };
    AuditoryPipeline pipeline;
    EXPECT_FALSE(pipeline.Init(models));
  }
}

TEST(AuditoryPipelineTest, CheckInitFromList) {
  AuditoryPipeline pipeline;
  const std::vector<AuditoryModelType> models = {
    MODEL_GAMMATONE_SLANEY, MODEL_MEDDIS_SYNAPSE_1986, MODEL_JACKSON_SPIKES };
  EXPECT_TRUE(pipeline.Init(models));
  EXPECT_EQ(models.size(), pipeline.Size());
}

TEST(AuditoryPipelineTest, CheckInitFromConfig) {
  AuditoryPipeline pipeline;
  AuditoryPipelineConfig config;
  config.add_models()->set_model_type(MODEL_GAMMATONE_SLANEY);
  config.add_models()->set_model_type(MODEL_MEDDIS_SYNAPSE_1986);
  config.add_models()->set_model_type(MODEL_JACKSON_SPIKES);
  EXPECT_TRUE(pipeline.Init(config));
  EXPECT_EQ(config.models().size(), pipeline.Size());
}

TEST(AuditoryPipelineTest, CheckStimulusFeed) {
  // Create pipeline config.
  AuditoryPipelineConfig config;
  config.add_models()->set_model_type(MODEL_GAMMATONE_SLANEY);
  config.add_models()->set_model_type(MODEL_MEDDIS_SYNAPSE_1986);

  // Initialize pipeline.
  AuditoryPipeline pipeline;
  EXPECT_TRUE(pipeline.Init(config));

  // Prepare stimulus.
  StimulusConfig stimulus_config;
  stimulus_config.set_sample_rate(16000 /* Hz */);
  constexpr int kNumEars = 1;
  constexpr int kNumSamples = 20000;
  const Eigen::ArrayXXd input = Eigen::ArrayXXd::Random(kNumEars, kNumSamples);

  // Run the stimulus through the pipeline. The output corresponds to the
  // synaptic firing rates.
  AuditoryResponse response;
  pipeline.Process(stimulus_config, input, &response);
  EXPECT_EQ(1, response.outputs().size());
  const auto &last_output = response.outputs()[0];
  EXPECT_EQ(OUTPUT_SYNAPSE_FIRING_RATES, last_output.first);
  const Eigen::ArrayXXd &output = last_output.second;
  EXPECT_EQ(kNumSamples, output.cols());
  EXPECT_LT(10, output.rows());  // Frequency channels.
  EXPECT_LT(5.0, output.minCoeff());  // Min firing rate.
  EXPECT_LT(50.0, output.maxCoeff());  // Max firing rate.
  const ChannelProperties &channel_properties = response.channel_properties();
  const auto &center_frequencies = channel_properties.center_frequencies;
  EXPECT_LT(0, center_frequencies.size());
  EXPECT_TRUE(std::is_sorted(center_frequencies.begin(),
                             center_frequencies.end()));

  // Run the same stimulus again, but this time store the intermediate outputs.
  stimulus_config.set_store_intermediate_outputs(true);
  pipeline.Process(stimulus_config, input, &response);
  const auto &outputs = response.outputs();
  EXPECT_EQ(2, outputs.size());
  EXPECT_EQ(OUTPUT_BASILAR_MEMBRANE_DISPLACEMENT, outputs[0].first);
  EXPECT_LT(0, outputs[0].second.rows());
  EXPECT_LT(0, outputs[0].second.cols());
  EXPECT_EQ(OUTPUT_SYNAPSE_FIRING_RATES, outputs[1].first);
  EXPECT_LT(0, outputs[1].second.rows());
  EXPECT_LT(0, outputs[1].second.cols());
}

}  // namespace audition
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
