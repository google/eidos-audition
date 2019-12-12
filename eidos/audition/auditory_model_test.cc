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

// Sanity checks for the auditory model interface.

#include "Eigen/Core"
#include "eidos/audition/auditory_model.h"
#include "eidos/audition/auditory_model_config.pb.h"
#include "eidos/audition/stimulus_config.pb.h"
#include "gtest/gtest.h"

namespace eidos {
namespace audition {
namespace {

// Model type (ID).
constexpr AuditoryModelType kModelType = MODEL_CARFAC;

// Auditory stage.
constexpr AuditoryStageType kStageType = STAGE_AUDITORY_NERVE_SYNAPSE;

// Model output.
constexpr AuditoryOutputType kOutputType = OUTPUT_SYNAPSE_FIRING_RATES;

// Mock a no-op synapse model.
class MockSynapseModel : public AuditoryModel {
 public:
  MockSynapseModel() : AuditoryModel(kModelType, kStageType, kOutputType,
                                     true /* waveform input */) {}

  ChannelProperties GetChannelProperties() const override { return {}; }

  void ProcessSegment(const Eigen::Ref<const Eigen::ArrayXXd> &input,
                      Eigen::ArrayXXd *output) override {}
};

}  // namespace

TEST(AuditoryModelTest, CheckMockModel) {
  MockSynapseModel model;
  EXPECT_EQ(kStageType, model.stage_type());
  EXPECT_EQ(kOutputType, model.output_type());
  EXPECT_EQ(kStageType, OutputType2Stage(kOutputType));
  EXPECT_TRUE(model.previous_model() == nullptr);
  EXPECT_EQ("MODEL_CARFAC", model.ModelTypeAsString());
}

}  // namespace audition
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
