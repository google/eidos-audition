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

// Basic waveform I/O check.

#include <vector>

#include "absl/strings/str_cat.h"
#include "gtest/gtest.h"
#include "wave/error.h"
#include "wave/file.h"

namespace eidos {

// Sample waveform.
const char kWavePath[] = "build/testdata/Untitled3.wav";

// Sample rate.
constexpr int kSampleRate = 44100;

TEST(WaveTest, Trivial) {
  wave::File file;
  ASSERT_EQ(wave::kNoError, file.Open(kWavePath, wave::OpenMode::kIn));
  EXPECT_EQ(file.sample_rate(), kSampleRate);
  EXPECT_EQ(file.bits_per_sample(), 16);
  EXPECT_EQ(file.channel_number(), 2);

  std::vector<float> contents;
  EXPECT_EQ(wave::kNoError, file.Read(&contents));
  EXPECT_EQ(contents.size() / file.channel_number(),
            5.558344671201814 * kSampleRate);
}

}  // namespace eidos

// Local Variables:
// mode: c++
// End:
