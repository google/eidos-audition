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

// Simple check that the auditory output wrapper behaves.

#include "eidos/audition/auditory_response.h"

#include <filesystem>
#include <cstdio>
#include <string>

#include "Eigen/Core"
#include "absl/strings/str_cat.h"
#include "eidos/audition/auditory_model_config.pb.h"
#include "gtest/gtest.h"

namespace eidos {
namespace audition {

constexpr int kNumChannels = 200;
constexpr int kNumSamples = 10000;

TEST(AuditoryResponseTest, CheckWriteNpy) {
  AuditoryResponse response;
  response.mutable_outputs()->emplace_back(
      OUTPUT_SYNAPSE_SPIKE_TIMES,
      Eigen::ArrayXXd::Random(kNumChannels, kNumSamples));
  const std::string temp_path = std::tmpnam(nullptr);
  response.WriteNpy(temp_path);
  EXPECT_TRUE(std::filesystem::remove(temp_path));
}

TEST(AuditoryResponseTest, CheckWriteNpz) {
  AuditoryResponse response;
  const Eigen::ArrayXXd &data = Eigen::ArrayXXd::Random(
      kNumChannels, kNumSamples);
  response.mutable_outputs()->emplace_back(
      OUTPUT_BASILAR_MEMBRANE_DISPLACEMENT, data);
  response.mutable_outputs()->emplace_back(
      OUTPUT_IHC_TRANSMEMBRANE_POTENTIAL, data);
  response.mutable_outputs()->emplace_back(
      OUTPUT_SYNAPSE_FIRING_RATES, data);
  response.mutable_outputs()->emplace_back(
      OUTPUT_SYNAPSE_SPIKE_TIMES, data);
  response.mutable_channel_properties()->center_frequencies.
      emplace_back(0.0);
  const std::string temp_path = std::tmpnam(nullptr);
  response.WriteNpz(temp_path);
  EXPECT_TRUE(std::filesystem::remove(temp_path));
}

}  // namespace audition
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
