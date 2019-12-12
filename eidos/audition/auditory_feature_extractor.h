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

// Offline auditory feature extractor supporting various models of auditory
// periphery.

#ifndef EIDOS_AUDITION_AUDITORY_FEATURE_EXTRACTOR_H_
#define EIDOS_AUDITION_AUDITORY_FEATURE_EXTRACTOR_H_

#include "eidos/audition/auditory_feature_extractor_config.pb.h"
#include "eidos/audition/auditory_response.h"
#include "eidos/audition/waveform.pb.h"

namespace eidos {
namespace audition {

class AuditoryFeatureExtractor {
 public:
  AuditoryFeatureExtractor() = default;
  ~AuditoryFeatureExtractor() = default;

  // Extracts filterbank features from the supplied mono (single-channel)
  // waveform and stores them in a matrix of dimension N by M, where
  // N is the number of filterbank coefficients and M is the number
  // of samples.
  //
  // The features may represent Neural Activity Patterns (NAP), IHC excitation,
  // Basilar Membrane (BM) displacements and so on.
  //
  // Given a filterbank that is producing outputs on a per-sample
  // basis, we want to downsample to produce outputs on a per-frame
  // basis.  In other words, given that the filterbank has produced
  // outputs for all the samples in the frame, we would like to
  // produce a single output [for each channel in the filterbank] that
  // summarizes the frame.
  //
  // Returns false if processing failed, true otherwise.
  bool ComputeFeatures(const Waveform &wave_data,
                       const AuditoryFeatureExtractorConfig &config,
                       AuditoryResponse *response) const;
};

}  // namespace audition
}  // namespace eidos

#endif  // EIDOS_AUDITION_AUDITORY_FEATURE_EXTRACTOR_H_

// Local Variables:
// mode: c++
// End:
