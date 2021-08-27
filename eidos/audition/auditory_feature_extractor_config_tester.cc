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

// Simple command-line utility for testing the supplied auditory feature
// extractor configuration.

#include <string>

#include "Eigen/Core"
#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "eidos/audition/auditory_feature_extractor.h"
#include "eidos/audition/auditory_feature_extractor_config.pb.h"
#include "eidos/audition/auditory_response.h"
#include "eidos/audition/waveform.pb.h"
#include "eidos/port/status-matchers.h"
#include "eidos/utils/proto_utils.h"
#include "gtest/gtest.h"

ABSL_FLAG(std::string, config_proto_file, "",
          "Configuration protocol buffer file in text format.");

namespace eidos {
namespace audition {

TEST(AuditoryFeatureExtractorConfigTest, CheckRunConfig) {
  // Read configuration.
  const std::string config_proto_file = absl::GetFlag(FLAGS_config_proto_file);
  ASSERT_TRUE(!config_proto_file.empty())
      << "Specify --config_proto_file";
  AuditoryFeatureExtractorConfig config;
  ASSERT_OK(utils::ReadTextProto(config_proto_file, &config))
      << "Failed to read configuration from: " << config_proto_file;
  config.mutable_config()->set_downsample_step(1);  // No downsampling.

  // Prepare the stimulus. All of the stimulus parameters are hardcoded at
  // present. The amplitude of random stimulus is unit normalized.
  Waveform input;
  input.set_num_channels(1);
  input.set_sample_rate(100000);  // 100 kHz.
  input.set_bits_per_sample(16);  // PCM.
  constexpr int kNumSamples = 50000;
  const Eigen::VectorXd vec = Eigen::VectorXd::Random(kNumSamples);
  input.mutable_samples()->Reserve(kNumSamples);
  for (int i = 0; i < kNumSamples; ++i) {
    input.add_samples(vec(i));
  }

  // Process the stimulus and perform some basic checks.
  AuditoryFeatureExtractor processor;
  AuditoryResponse response;
  EXPECT_TRUE(processor.ComputeFeatures(input, config, &response))
      << "Failed to compute features";
  const auto &outputs = response.outputs();
  EXPECT_EQ(1, outputs.size());
  EXPECT_EQ(kNumSamples, outputs[0].second.cols());
}

}  // namespace audition
}  // namespace eidos

int main(int argc, char *argv[]) {
  absl::ParseCommandLine(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

// Local Variables:
// mode: c++
// End:
