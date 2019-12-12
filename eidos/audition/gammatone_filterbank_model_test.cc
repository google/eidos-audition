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

// Simple test for gammatone filterbanks.

#include <algorithm>
#include <functional>

#include "Eigen/Core"
#include "eidos/audition/stimulus_config.pb.h"
#include "gtest/gtest.h"
#include "third_party/audition/models/gammatone/gammatone_filterbank_model.h"

namespace eidos {
namespace audition {
namespace {

// A small amplitude value.
constexpr float kSmallAmplitude = 1E-6;

// Checks computed channel properties. The Slaney's implementation places center
// frequencies in decreasing order, while Ma's implementation uses increasing
// order.
template <typename Comparator>
void CheckChannelProperties(const ChannelProperties &channel_properties,
                            int expected_num_channels,
                            Comparator comparator) {
  const auto &center_frequencies = channel_properties.center_frequencies;
  EXPECT_EQ(center_frequencies.size(), expected_num_channels);
  EXPECT_TRUE(std::is_sorted(center_frequencies.begin(),
                             center_frequencies.end(), comparator));
}

// Checks random unit impulse response of Slaney's filterbank.
void CheckUnitImpulseResponseSlaney(int num_channels, int num_expected_channels,
                                    float min_cf_hz, float expected_min_cf_hz,
                                    float max_cf_hz, float expected_max_cf_hz) {
  constexpr int kNumEars = 1;
  constexpr int kNumSamples = 20000;  // Long enough so we get a good average.
  Eigen::ArrayXXd input = Eigen::ArrayXXd::Zero(kNumEars, kNumSamples);
  input(0, 1) = 1.0;  // Unit impulse.

  // Set up stimulus configuration.
  constexpr int kSampleRate = 16000;  // Hz.
  StimulusConfig stimulus_config;
  stimulus_config.set_sample_rate(kSampleRate);
  stimulus_config.set_num_channels(num_channels);
  stimulus_config.set_lowest_cf_hz(min_cf_hz);
  stimulus_config.set_highest_cf_hz(max_cf_hz);

  // Initialize the model and process random stimulus.
  GammatoneFilterbankConfig config;
  config.set_filter_type(GAMMATONE_FILTERBANK_SLANEY);
  GammatoneFilterbankModel model(config);
  model.Init(stimulus_config);
  Eigen::ArrayXXd output;
  model.ProcessSegment(input, &output);

  // Check channel properties. Note: The frequencies extend in descending order.
  // The lowest and highest CFs of the filters don't exactly correspond to the
  // requested values due to specific ERB computations in Slaney's filter.
  const ChannelProperties &channel_properties = model.GetChannelProperties();
  const auto &cfs = channel_properties.center_frequencies;
  CheckChannelProperties(channel_properties, output.rows(),
                         std::greater<double>());
  EXPECT_NEAR(expected_min_cf_hz, cfs.back(), 1E-6);
  EXPECT_LT(2, cfs.size());
  const double cf_delta = cfs[cfs.size() - 2] - cfs.back();
  EXPECT_LT(0.0, cf_delta);
  EXPECT_LT(cfs[0] + cf_delta / 2.0, expected_max_cf_hz);

  // Check response properties.
  EXPECT_EQ(num_expected_channels, output.rows());
  EXPECT_EQ(kNumSamples, output.cols());
  const int total_samples = channel_properties.NumChannels() * kNumSamples;
  EXPECT_EQ(output.size(), total_samples);
  EXPECT_LT(-1.0, output.minCoeff());
  EXPECT_GT(0.5, output.maxCoeff());
  const float output_amplification_level = output.matrix().squaredNorm();
  EXPECT_LT(0.8, output_amplification_level);

  // Compute percentage of the amplitudes we can safely neglect.
  int negligible_samples = 0;
  for (int i = 0; i < output.rows(); ++i) {  // channels
    for (int j = 0; j < kNumSamples; ++j) {
      if (output(i, j) < kSmallAmplitude) negligible_samples++;
    }
  }
  const float negligible_samples_ratio =
      negligible_samples / static_cast<float>(total_samples);
  EXPECT_LT(0.95, negligible_samples_ratio);
}

}  // namespace

TEST(GammatoneFilterbankModelTest, UnitImpulseResponseSlaneyDefNumChannels) {
  CheckUnitImpulseResponseSlaney(
      0, GammatoneFilterbankModel::kDefaultNumChannels,  // Number of channels.
      0, GammatoneFilterbankModel::kLowestCFrequency,  // Lowest CF.
      // Expect maximum CF to be around half the sample rate.
      0, 8000.0);  // Highest CF.
}

TEST(GammatoneFilterbankModelTest, CheckUnitImpulseResponseSlaney) {
  constexpr int kRequestedChannels = 155;
  constexpr double kMinCfHz = 100.0;
  constexpr double kMaxCfHz = 5000.0;
  CheckUnitImpulseResponseSlaney(
      kRequestedChannels, kRequestedChannels,
      kMinCfHz, kMinCfHz,
      kMaxCfHz, kMaxCfHz);
}

namespace {

// Checks gammatone response given the implementation on a random stimulus.
template <typename Comparator>
void CheckRandomResponse(GammatoneFilterbankType filter_type,
                         float min_displacement,
                         float max_displacement,
                         Comparator comparator) {
  constexpr double kSampleRate = 44000;  // Hz.
  StimulusConfig stimulus_config;
  stimulus_config.set_sample_rate(kSampleRate);

  constexpr int kNumEars = 1;
  constexpr int kNumSamples = 50000;  // Long enough so we get a good average.
  const Eigen::ArrayXXd input = Eigen::ArrayXXd::Random(kNumEars, kNumSamples);

  GammatoneFilterbankConfig config;
  config.set_filter_type(filter_type);
  GammatoneFilterbankModel model(config);
  model.Init(stimulus_config);
  Eigen::ArrayXXd output;
  model.ProcessSegment(input, &output);
  const ChannelProperties &channel_properties = model.GetChannelProperties();
  CheckChannelProperties(channel_properties, output.rows(), comparator);
  EXPECT_EQ(kNumSamples, output.cols());
  const int total_samples = channel_properties.NumChannels() * kNumSamples;
  EXPECT_EQ(output.size(), total_samples);
  EXPECT_LT(min_displacement, output.minCoeff());
  EXPECT_GT(max_displacement, output.maxCoeff());
  const float output_amplification_level = output.matrix().squaredNorm();
  EXPECT_LT(0.8, output_amplification_level);
}

}  // namespace

TEST(GammatoneFilterbankModelTest, CheckRandomResponseSlaney) {
  CheckRandomResponse(GAMMATONE_FILTERBANK_SLANEY,
                      -1.0,  // min_displacement
                      0.8,  // max_displacement
                      std::greater<double>());
}

TEST(GammatoneFilterbankModelTest, CheckRandomResponseCookeAndMa) {
  CheckRandomResponse(GAMMATONE_FILTERBANK_COOKE_AND_MA,
                      -520.0,  // min_displacement
                      520.0,  // max displacement
                      std::less<double>());
}

}  // namespace audition
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
