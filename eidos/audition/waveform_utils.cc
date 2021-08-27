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

#include "eidos/audition/waveform_utils.h"

#include <cstdint>
#include <limits>
#include <vector>

#include "eidos/port/logging.h"
#include "google/protobuf/repeated_field.h"
#include "wave/error.h"
#include "wave/file.h"

namespace eidos {
namespace audition {
namespace {

// Scaling factor.
const float kInt16ToFloatDivisor =
    -static_cast<float>(std::numeric_limits<int16_t>::min());

// Converts wave I/O error codes to readable strings.
std::string ErrorToString(wave::Error status) {
  switch (status) {
    case wave::kFailedToOpen:
      return "Open failed";
    case wave::kNotOpen:
      return "File not open";
    case wave::kInvalidFormat:
      return "Invalid format";
    case wave::kWriteError:
      return "Write error";
    case wave::kReadError:
      return "Read error";
    case wave::kInvalidSeek:
      return "Invalid seek";
    default:
      return "Unknown";
  }
}

}  // namespace

bool ReadWaveform(const std::string &filepath, Waveform *proto) {
  wave::File file;
  wave::Error status = file.Open(filepath, wave::OpenMode::kIn);
  if (status != wave::kNoError) {
    GOOGLE_LOG(ERROR) << "Failed to open waveform from \"" << filepath << "\": "
                      << ErrorToString(status);
    return false;
  }
  proto->Clear();
  proto->set_num_channels(file.channel_number());
  if (proto->num_channels() != 1) {
    GOOGLE_LOG(ERROR) << "Multi-channel audio is not supported yet";
    return false;
  }
  proto->set_sample_rate(file.sample_rate());
  proto->set_bits_per_sample(file.bits_per_sample());

  std::vector<float> contents;
  status = file.Read(&contents);
  if (status != wave::kNoError) {
    GOOGLE_LOG(ERROR) << "Failed to read waveform from \"" << filepath << "\": "
                      << ErrorToString(status);
  }
  *proto->mutable_samples() = {contents.begin(), contents.end()};
  if (proto->samples().size() == 0) {
    GOOGLE_LOG(ERROR) << "Waveform is empty: " << filepath;
    return false;
  }
  GOOGLE_LOG(INFO) << "Waveform info: Sample rate: " << proto->sample_rate()
                   << " Hz, samples: " << proto->samples().size()
                   << ", channels: " << proto->num_channels()
                   << ", bits per sample: " << proto->bits_per_sample();
  return true;
}

bool WriteWaveform(const Waveform &proto, const std::string &filepath) {
  wave::File file;
  wave::Error status = file.Open(filepath, wave::OpenMode::kOut);
  if (status != wave::kNoError) {
    GOOGLE_LOG(ERROR) << "Failed to open waveform in " << filepath << ": "
                      << ErrorToString(status);
    return false;
  }
  file.set_channel_number(proto.num_channels());
  file.set_sample_rate(proto.sample_rate());
  file.set_bits_per_sample(proto.bits_per_sample());
  status = file.Write({proto.samples().begin(), proto.samples().end()});
  if (status != wave::kNoError) {
    GOOGLE_LOG(ERROR) << "Failed to write waveform to " << filepath << ": "
                      << ErrorToString(status);
  }
  return true;
}

void ScaleToUnitWaveform(const Waveform &waveform, float scaling_gain,
                         std::vector<float> *samples) {
  const int num_samples = waveform.samples().size();
  samples->clear();
  samples->resize(num_samples);
  for (int i = 0; i < num_samples; ++i) {
    (*samples)[i] = (waveform.samples(i) / kInt16ToFloatDivisor) * scaling_gain;
  }
}

namespace {

template <typename T>
void DenormalizeWaveformImpl(Waveform *waveform) {
  std::transform(
      waveform->mutable_samples()->begin(), waveform->mutable_samples()->end(),
      waveform->mutable_samples()->begin(), [](float sample) -> float {
        return sample * std::numeric_limits<T>::max();
      });
}

}  // namespace

void DenormalizeWaveform(Waveform *waveform) {
  GOOGLE_CHECK(waveform != nullptr);
  switch (waveform->bits_per_sample()) {
    case 8:
      DenormalizeWaveformImpl<int8_t>(waveform);
      break;
    case 16:
      DenormalizeWaveformImpl<int16_t>(waveform);
      break;
    case 32:
      DenormalizeWaveformImpl<int32_t>(waveform);
      break;
    default: {
      GOOGLE_LOG(FATAL) << "Unsupported sample width: "
                        << waveform->bits_per_sample();
    }
  }
}

}  // namespace audition
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
