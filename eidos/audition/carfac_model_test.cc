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

// Simple test for CARFAC model.

#include <algorithm>
#include <cstddef>
#include <functional>

#include "Eigen/Core"
#include "eidos/audition/carfac_config.pb.h"
#include "eidos/audition/carfac_model.h"
#include "eidos/audition/stimulus_config.pb.h"
#include "eidos/stubs/logging.h"
#include "gtest/gtest.h"

namespace eidos {
namespace audition {
namespace {

// Given the model configuration runs the filterbank over the simple stimulus
// and checks the resulting signal. The expected level of the output signal is
// provided by <min_aplification_level>.
void BasicFilterbankCheck(const CarfacConfig &model_config,
                          int min_amplification_level) {
  // Generate simple random stimulus.
  StimulusConfig stimulus_config;
  stimulus_config.set_audio_scaling_gain(1.0);
  stimulus_config.set_sample_rate(16000 /* 16 kHz */);

  const int kNumEars = 1;
  const int kNumSamples = 10000;  // Long enough so we get a good average.
  const Eigen::ArrayXXd input = Eigen::ArrayXXd::Random(kNumEars, kNumSamples);

  // Filter the stimulus.
  CarfacModel model(model_config);
  model.Init(stimulus_config);
  Eigen::ArrayXXd output;
  model.ProcessSegment(input, &output);

  // Perform basic sanity checks.
  const ChannelProperties &channel_properties = model.GetChannelProperties();
  EXPECT_EQ(kNumSamples, output.cols());
  EXPECT_LT(0, channel_properties.NumChannels());
  EXPECT_EQ(output.rows(), channel_properties.NumChannels());
  const auto &center_frequencies = channel_properties.center_frequencies;
  // Note: center frequencies are sorted in descending order.
  EXPECT_TRUE(std::is_sorted(center_frequencies.begin(),
                             center_frequencies.end(), std::greater<double>()));
  const float output_level = output.matrix().squaredNorm();
  GOOGLE_LOG(INFO) << "Output level: " << output_level;
  EXPECT_LT(min_amplification_level, output_level);
}

}  // namespace

// Performs generation of neural activity pattern (NAP) triggered by default
// configuration, accompanied by some basic checks.
TEST(CarfacModelTest, SimpleNapTest) {
  CarfacConfig config;
  constexpr double kMinAmplificationLevel = 1E+5;
  BasicFilterbankCheck(config, kMinAmplificationLevel);
}

// Obtains estimates of basilar membrane (BM) displacements and performs simple
// checks.
TEST(CarfacModelTest, SimpleBmTest) {
  CarfacConfig config;
  config.mutable_output_types()->set_store_nap(false);
  config.mutable_output_types()->set_store_bm(true);
  constexpr double kMinAmplificationLevel = 1E+5;
  BasicFilterbankCheck(config, kMinAmplificationLevel);
}

// Obtains outer hair cell (OHC) signal and performs simple checks.
TEST(CarfacModelTest, SimpleOhcTest) {
  CarfacConfig config;
  config.mutable_output_types()->set_store_nap(false);
  config.mutable_output_types()->set_store_ohc(true);
  constexpr double kMinAmplificationLevel = 1E+5;
  BasicFilterbankCheck(config, kMinAmplificationLevel);
}

// Obtains adaptive gain control (AGC) signal and performs simple checks.
TEST(CarfacModelTest, SimpleAgcTest) {
  CarfacConfig config;
  config.mutable_output_types()->set_store_nap(false);
  config.mutable_output_types()->set_store_agc(true);
  constexpr double kMinAmplificationLevel = 2E+3;
  BasicFilterbankCheck(config, kMinAmplificationLevel);
}

// Disables AGC through the coniguration.
TEST(CarfacModelTest, SimpleNapNoAgcTest) {
  CarfacConfig config;
  config.mutable_agc()->set_num_stages(0);  // Disable AGC.
  constexpr double kMinAmplificationLevel = 1E+5;
  BasicFilterbankCheck(config, kMinAmplificationLevel);
}

// Tweaks some of the Cascade of Asymmetric Resonators (CAR) filter parameters
// and checks that the model is still sane.
TEST(CarfacModelTest, SimpleTweakCarTest) {
  CarfacConfig config;
  constexpr double kMinAmplificationLevel = 1E+5;
  config.mutable_car()->set_min_pole_hz(20.0);  // Hz.
  for (int i = 1; i <= 10; ++i) {
    config.mutable_car()->set_erb_per_step(0.1);
    BasicFilterbankCheck(config, kMinAmplificationLevel);
  }
}

// CARFAC's number of analysis channels depends on the sampling rate.
TEST(CarfacModelTest, SamplingRatesTest) {
  const int kNumEars = 1;
  const int kNumSamples = 10000;  // Long enough so we get a good average.
  const Eigen::ArrayXXd input = Eigen::ArrayXXd::Random(kNumEars, kNumSamples);
  const std::vector<float> sampling_rates_hz = {
    8000, 16000, 22040, 44100, 80000, 100000};
  const std::vector<int> num_channels = { 53, 65, 71, 83, 94, 98 };
  for (std::size_t i = 0; i < sampling_rates_hz.size(); ++i) {
    StimulusConfig stimulus_config;
    const int sample_rate = sampling_rates_hz[i];
    stimulus_config.set_sample_rate(sample_rate);

    CarfacConfig model_config;
    CarfacModel model(model_config);
    model.Init(stimulus_config);
    Eigen::ArrayXXd output;
    model.ProcessSegment(input, &output);
    EXPECT_EQ(output.rows(), num_channels[i]);
  }
}

}  // namespace audition
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
