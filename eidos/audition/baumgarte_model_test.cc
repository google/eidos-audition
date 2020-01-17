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

// Simple test for Baumgarte's model of auditory periphery.

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <functional>
#include <string>

#include "Eigen/Core"
#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"
#include "eidos/audition/stimulus_config.pb.h"
#include "eidos/audition/waveform.pb.h"
#include "eidos/audition/waveform_utils.h"
#include "eidos/stubs/logging.h"
#include "gtest/gtest.h"
#include "third_party/audition/models/baumgarte/baumgarte_model.h"

namespace eidos {
namespace audition {
namespace {

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

}  // namespace

TEST(BaumgarteModelTest, SimpleTest) {
  StimulusConfig stimulus_config;
  stimulus_config.set_sample_rate(100000 /* Hz */);
  stimulus_config.set_audio_scaling_gain(1.0);
  constexpr int kDownsampleStep = 5;
  stimulus_config.set_downsample_step(kDownsampleStep);

  constexpr int kNumEars = 1;
  constexpr int kNumSamples = 10000;  // Long enough so we get a good average.
  const Eigen::ArrayXXd input = Eigen::ArrayXXd::Random(kNumEars, kNumSamples);

  BaumgarteModel model;
  model.Init(stimulus_config);
  Eigen::ArrayXXd output;
  model.ProcessSegment(input, &output);
  EXPECT_EQ(output.cols(), kNumSamples / kDownsampleStep);
  const ChannelProperties &channel_properties = model.GetChannelProperties();
  CheckChannelProperties(channel_properties, output.rows(), std::less<double>());
  EXPECT_EQ(output.size(), channel_properties.NumChannels() *
            kNumSamples / kDownsampleStep);
  const float output_level = output.matrix().squaredNorm();
  GOOGLE_LOG(INFO) << "Output level: " << output_level;
  constexpr double kMinOutputAmplificationLevel = 2.0;
  EXPECT_LT(kMinOutputAmplificationLevel, output_level);
  EXPECT_LT(-1.0, output.minCoeff());  // Minimal potential.
  EXPECT_GT(2.0, output.maxCoeff());   // Maximal potential.
}

// Checks various values for downsampling count parameter.
TEST(BaumgarteModelTest, DownsampleTest) {
  constexpr int kNumEars = 1;
  constexpr int kNumSamples = 10000;
  const Eigen::ArrayXXd input = Eigen::ArrayXXd::Random(kNumEars, kNumSamples);

  constexpr int kMaxDownsample = 20;
  constexpr int kDownsampleStep = 4;
  StimulusConfig stimulus_config;
  stimulus_config.set_sample_rate(100000 /* Hz */);
  stimulus_config.set_audio_scaling_gain(1.0);
  for (int i = 1; i < kMaxDownsample; i += kDownsampleStep) {
    stimulus_config.set_downsample_step(i);
    BaumgarteModel model;
    model.Init(stimulus_config);
    Eigen::ArrayXXd output;
    model.ProcessSegment(input, &output);
    const int expected_num_samples = std::ceil(
        static_cast<double>(kNumSamples) / i);
    EXPECT_EQ(output.cols(), expected_num_samples);
  }
}

// Checks various values for the number of analysis channels.
TEST(BaumgarteModelTest, NumberOfChannelsTest) {
  constexpr int kNumEars = 1;
  constexpr int kNumSamples = 10000;
  const Eigen::ArrayXXd input = Eigen::ArrayXXd::Random(kNumEars, kNumSamples);

  constexpr int kMaxNumChannels = 280;
  constexpr int kMinNumChannels = 100;
  constexpr int kStep = 5;
  StimulusConfig stimulus_config;
  stimulus_config.set_sample_rate(100000 /* Hz */);
  stimulus_config.set_audio_scaling_gain(1.0);
  for (int i = kMinNumChannels; i <= kMaxNumChannels; i += kStep) {
    stimulus_config.set_num_channels(i);
    BaumgarteModel model;
    model.Init(stimulus_config);
    Eigen::ArrayXXd output;
    model.ProcessSegment(input, &output);
    EXPECT_EQ(output.rows(), i);
  }
}

namespace {

// Path the test data for auditory filters.
const char kAuditoryTestDataPath[] = "eidos/audition/testdata/auditory";

// Number of reference Basilar Membrane (BM) sections.
constexpr int kRefNumBmSections = 251;

// Number of reference time steps.
constexpr int kRefNumTimeSteps = 1000;

// Returns full path pointing to the auditory test file.
std::string AuditoryTestFile(absl::string_view filename) {
  return absl::StrCat(kAuditoryTestDataPath, "/", filename);
}

// Utility function to read a vector of specified POD types from a headerless
// binary file.
template <typename T>
void ReadFromFile(const std::string &path, std::vector<T> *vec) {
  std::ifstream ifs(path, std::ios::binary);
  ASSERT_TRUE(ifs.good()) << "Failed to open " << path;
  ifs.seekg(0, std::ios::end);
  const std::size_t file_size = ifs.tellg();
  ifs.seekg(0, std::ios::beg);
  const int num_values = file_size / sizeof(T);
  EXPECT_EQ(kRefNumBmSections * kRefNumTimeSteps, num_values);
  vec->resize(num_values);
  ifs.read(reinterpret_cast<char *>(vec->data()), file_size);
  EXPECT_TRUE(ifs.good()) << "Failed to read " << num_values
                          << " float values from " << path;
}

// Runs the stimulus waveform through the filterbank and compares the resulting
// excitation against the golden data from a file.
void CheckExcitationImpl(const std::string &waveform_path,
                         const std::string &ref_excitation_path) {
  // Read the waveform.
  Waveform waveform;
  ASSERT_TRUE(ReadWaveform(waveform_path, &waveform))
      << "Failed to open " << waveform_path;
  DenormalizeWaveform(&waveform);
  const int kNumEars = 1;  // Mono.
  const int num_input_samples = waveform.samples().size();
  Eigen::ArrayXXd input(kNumEars, num_input_samples);
  for (int i = 0; i < num_input_samples; ++i) {
    input(0, i) = waveform.samples(i);
  }

  // Generate excitation. Use downsample step of 10.
  StimulusConfig stimulus_config;
  stimulus_config.set_audio_scaling_gain(1.0);
  stimulus_config.set_sample_rate(waveform.sample_rate());
  stimulus_config.set_downsample_step(10);
  BaumgarteModel model;
  model.Init(stimulus_config);
  Eigen::ArrayXXd output;
  model.ProcessSegment(input, &output);
  const ChannelProperties &channel_properties = model.GetChannelProperties();
  CheckChannelProperties(channel_properties, output.rows(), std::less<double>());
  EXPECT_EQ(kRefNumTimeSteps, output.cols());
  const Eigen::ArrayXXd &exc = output;
  EXPECT_EQ(kRefNumBmSections, exc.rows());
  EXPECT_EQ(kRefNumTimeSteps, exc.cols());

  // Read the reference excitation. These are stored in an N-by-M matrix, where
  // the N corresponds to BM sections and M corresponds to time.
  std::vector<float> ref_exc_vals;
  ReadFromFile(ref_excitation_path, &ref_exc_vals);
  const Eigen::Map<Eigen::Matrix<float, kRefNumBmSections, kRefNumTimeSteps>>
      ref_exc_map(ref_exc_vals.data());
  const Eigen::ArrayXXd ref_exc = ref_exc_map.cast<double>();
  EXPECT_EQ(kRefNumBmSections, ref_exc.rows());
  EXPECT_EQ(kRefNumTimeSteps, ref_exc.cols());

  // We almost match the reference implementation.
  const float kEpsilon = 7.9E-3;
  EXPECT_TRUE(ref_exc.isApprox(exc, kEpsilon));
}

// Runs the stimulus waveform through the filterbank and compares the resulting
// excitation against the golden data from a file.
void CheckExcitation(const std::string &stimulus_filename,
                     const std::string &ref_excitation_filename) {
  CheckExcitationImpl(AuditoryTestFile(stimulus_filename),
                      AuditoryTestFile(ref_excitation_filename));
}

}  // namespace

// Checks response to pure 1kHz tone at sound pressure level (SPL) of 60 dB.
TEST(BaumgarteModelTest, Check_1kHz60dB_SPL_Tone) {
  CheckExcitation("1kHz_60dB_SPL_tone.wav", "1kHz_60dB_SPL_tone_baumgarte.exc");
}

// Checks response to two pure 1kHz impulses at 0s and 50ms.
TEST(BaumgarteModelTest, Check_PureImpulse) {
  CheckExcitation("impulse_1000at0s_and_-1000at0.05s.wav",
                  "impulse_1000at0s_and_-1000at0.05s_baumgarte.exc");
}

}  // namespace audition
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
