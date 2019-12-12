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

// Utilities for handling waveform files.

#ifndef EIDOS_AUDITION_WAVEFORM_UTIL_H_
#define EIDOS_AUDITION_WAVEFORM_UTIL_H_

#include <string>
#include <vector>

#include "eidos/audition/waveform.pb.h"

namespace eidos {
namespace audition {

// Given the waveform file path reads it into a proto. Returns false in case
// of error. The waveform samples are stored normalized in memory. In order
// to retrieve the original values please use "DenormalizeWaveform" API below.
bool ReadWaveform(const std::string &filepath, Waveform *proto);

// Saves the waveform to the specified file path. Returns false in case of
// error. The original waveform is presumed to be normalized in memory.
// Denormalization happens before saving it to disk.
bool WriteWaveform(const Waveform &proto, const std::string &filepath);

// Takes the waveform and scales it to [-1.0,1.0] range applying gain
// factor on the result.
void ScaleToUnitWaveform(const Waveform &waveform, float scaling_gain,
                         std::vector<float> *samples);

// By default the input waveform is normalized by the reader. Denormalize
// the waveform restoring the original values.
void DenormalizeWaveform(Waveform *waveform);

}  // namespace audition
}  // namespace eidos

#endif  // EIDOS_AUDITION_WAVEFORM_UTIL_H_

// Local Variables:
// mode: c++
// End:
