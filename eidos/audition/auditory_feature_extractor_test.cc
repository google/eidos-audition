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

// Sanity check for the feature extractor.

#include "Eigen/Core"
#include "eidos/audition/auditory_feature_extractor.h"
#include "eidos/audition/auditory_feature_extractor_config.pb.h"
#include "eidos/audition/auditory_response.h"
#include "eidos/audition/waveform_utils.h"
#include "gtest/gtest.h"

namespace eidos {
namespace audition {
namespace {

// Source waveform in original sample rate (24 kHz).
const char kWaveFilePath[] =
    "eidos/audition/testdata/ne_np/nep_letter_ma.wav";

// Source waveform in 100kHz.
const char k100kHzWaveFilePath[] =
    "eidos/audition/testdata/ne_np/nep_letter_ma_100kHz.wav";

// CARFAC filterbank produces 72 tonotopic channels for 24 kHz input.
const int kNumChannelsAt24Carfac = 72;

// Number of channels for Baumgarte's model.
const int kNumChannelsBaumgarte = 251;

}  // namespace

class AuditoryFeatureExtractorTest : public ::testing::Test {
 protected:
  // Some filterbanks (like Baumgarte's and Zilany's) require the stimulus
  // to be 100kHz or more.
  void Init(bool high_sample_rate) {
    if (high_sample_rate) {
      ASSERT_TRUE(ReadWaveform(k100kHzWaveFilePath, &wave_));
      EXPECT_EQ(100000, wave_.sample_rate());
    } else {
      ASSERT_TRUE(ReadWaveform(kWaveFilePath, &wave_));
      EXPECT_EQ(24000, wave_.sample_rate());
    }
    EXPECT_EQ(1, wave_.num_channels());  // Mono.
  }

  void CheckModel(AuditoryModelType model_type,
                  int expected_num_channels) {
    StimulusConfig *stimulus_config = config_.mutable_config();
    stimulus_config->set_audio_scaling_gain(0.07);
    config_.mutable_pipeline()->add_models()->set_model_type(model_type);
    AuditoryResponse response;
    AuditoryFeatureExtractor extractor;
    EXPECT_TRUE(extractor.ComputeFeatures(wave_, config_, &response))
        << "Failed to extract " << AuditoryModelType_Name(model_type)
        << " features!";
    EXPECT_EQ(1, response.outputs().size());
    const Eigen::ArrayXXd &output = response.outputs().back().second;
    EXPECT_EQ(wave_.samples().size(), output.cols());
    EXPECT_EQ(expected_num_channels, output.rows());
    // Check that output is not all zeros.
    EXPECT_GT(output.matrix().squaredNorm(), 0.0);
  }

  // Configuration for the cochlear extractor.
  AuditoryFeatureExtractorConfig config_;

  // Actual waveform for testing.
  Waveform wave_;
};

// Runs the input waveform through the filterbank only without any downsampling
// or audio image stabilization. The output is sample-based.
TEST_F(AuditoryFeatureExtractorTest, CarfacCheckOriginalWaveform) {
  Init(false /* use native sample rate */);
  CheckModel(MODEL_CARFAC, kNumChannelsAt24Carfac);
}

TEST_F(AuditoryFeatureExtractorTest, BaumgarteCheckOriginalWaveform) {
  Init(true /* use 100kHz audio */);
  config_.mutable_config()->set_downsample_step(1 /* no downsampling */);
  CheckModel(MODEL_BAUMGARTE, kNumChannelsBaumgarte);
}

TEST_F(AuditoryFeatureExtractorTest, Zilany2014CheckOriginalWaveform) {
  Init(true /* use 100kHz audio */);
  config_.mutable_config()->set_downsample_step(1 /* no downsampling */);
  CheckModel(
      MODEL_ZILANY_IHC_2014,
      kNumChannelsBaumgarte /* same number of channels as Baumgarte's */);
}

}  // namespace audition
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
