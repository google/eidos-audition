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

// Sanity check for waveform utilities.

#include <algorithm>
#include <cmath>
#include <numeric>
#include <string>
#include <vector>

#include "absl/strings/str_cat.h"
#include "eidos/audition/waveform_utils.h"
#include "eidos/utils/test_utils.h"
#include "google/protobuf/util/message_differencer.h"
#include "gtest/gtest.h"

using google::protobuf::util::MessageDifferencer;

namespace eidos {
namespace audition {
namespace {

// Source waveform.
const char kWaveFilePath[] =
    "eidos/audition/testdata/ne_np/nep_0258_0119737288.wav";

}  // namespace

TEST(WaveformUtilsTest, Reading) {
  Waveform wave;
  ASSERT_TRUE(ReadWaveform(kWaveFilePath, &wave));
  EXPECT_EQ(1, wave.num_channels());  // Mono.
  EXPECT_EQ(48000, wave.sample_rate());
  EXPECT_EQ(16, wave.bits_per_sample());
  EXPECT_EQ(192912, wave.samples().size());
}

TEST(WaveformUtilsTest, ReadingAndWriting) {
  Waveform wave;
  ASSERT_TRUE(ReadWaveform(kWaveFilePath, &wave));
  const auto temp_path_info = utils::TempPath();
  const std::string &temp_path = temp_path_info.second;
  ASSERT_TRUE(WriteWaveform(wave, temp_path));

  Waveform new_wave;
  ASSERT_TRUE(ReadWaveform(temp_path, &new_wave));
  EXPECT_TRUE(MessageDifferencer::Equals(new_wave, wave));
}

TEST(WaveformUtilsTest, UnitScaling) {
  Waveform wave;
  ASSERT_TRUE(ReadWaveform(kWaveFilePath, &wave));
  std::vector<float> samples;
  ScaleToUnitWaveform(wave, 0.5, &samples);
  EXPECT_LT(0, samples.size());
  for (auto sample : samples) {
    EXPECT_GT(1.0, sample);
    EXPECT_LT(-1.0, sample);
  }
}

TEST(WaveformUtilsTest, CheckDenormalize) {
  Waveform wave;
  ASSERT_TRUE(ReadWaveform(kWaveFilePath, &wave));
  for (auto sample : wave.samples()) {
    EXPECT_GT(1.0, std::fabs(sample));
  }
  DenormalizeWaveform(&wave);
  const float mean =
      std::accumulate(wave.samples().begin(), wave.samples().end(), 0.0) /
      wave.samples().size();
  EXPECT_LT(2.0, std::fabs(mean));  // Samples are not within [-1, 1] range.
}

}  // namespace audition
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
