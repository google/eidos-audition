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

#include "eidos/audition/auditory_feature_extractor.h"

#include <vector>

#include "eidos/audition/auditory_pipeline.h"
#include "eidos/audition/resample.h"
#include "eidos/audition/stimulus_config.pb.h"
#include "eidos/audition/waveform_utils.h"
#include "eidos/audition/windowing.h"
#include "eidos/stubs/logging.h"

namespace eidos {
namespace audition {
namespace {

// Only accept monaural input for now.
constexpr int kNumEars = 1;

}  // namespace

bool AuditoryFeatureExtractor::ComputeFeatures(
    const Waveform &waveform, const AuditoryFeatureExtractorConfig &config,
    AuditoryResponse *response) const {
  // Fetch the audio sampling rate from the waveform.
  StimulusConfig stimulus_config = config.config();
  stimulus_config.set_sample_rate(waveform.sample_rate());

  // Scales the input waveform to [-1, 1] range.
  std::vector<float> samples;
  ScaleToUnitWaveform(waveform, stimulus_config.audio_scaling_gain(), &samples);
  const int num_samples = samples.size();
  Eigen::ArrayXXd input(kNumEars, num_samples);
  for (int i = 0; i < num_samples; ++i) {
    input(0, i) = samples[i];
  }

  // Initialize the pipeline.
  AuditoryPipeline pipeline;
  if (!pipeline.Init(config.pipeline())) {
    GOOGLE_LOG(ERROR) << "Failed to initialize auditory pipeline";
    return false;
  }

  // Feed the stimulus to the pipeline.
  pipeline.Process(stimulus_config, input, response);
  const auto &last_output = response->outputs().back().second;
  if (last_output.matrix().squaredNorm() <= 0) {
    GOOGLE_LOG(ERROR) << "Invalid filterbank output!";
    return false;
  }
  if (stimulus_config.downsample_step() > 1) {
    // Adjust the sample rate by downsample step.
    stimulus_config.set_sample_rate(stimulus_config.sample_rate() /
                                    stimulus_config.downsample_step());
  }
  GOOGLE_LOG(INFO) << "Generated " << last_output.size() << " values "
                   << "[sample_rate: " << stimulus_config.sample_rate()
                   << ", frequency channels: " << last_output.rows()
                   << ", samples: " << last_output.cols() << "].";

  // Resample the output, if required.
  if (ShouldResampleOutputs(stimulus_config)) {
    for (auto &output : *response->mutable_outputs()) {
      Eigen::ArrayXXd output_frames = Resample(
          output.second, stimulus_config);
      output.second.swap(output_frames);
    }

    // Update sampling rate.
    const int output_sample_rate = stimulus_config.sample_rate() *
        stimulus_config.output_resample_up_factor() /
        stimulus_config.output_resample_down_factor();
    stimulus_config.set_sample_rate(output_sample_rate);
    GOOGLE_LOG(INFO) << "Resampler: Generated " << last_output.cols()
                     << " samples. Sample rate: "
                     << stimulus_config.sample_rate();
  }

  // Apply windowing if necessary.
  if (stimulus_config.apply_window_to_outputs()) {
    for (auto &output : *response->mutable_outputs()) {
      Eigen::ArrayXXd output_frames = WindowAndIntegrateTime(
          output.second, stimulus_config);
      output.second.swap(output_frames);
    }
    GOOGLE_LOG(INFO) << "Windowing: Generated " << last_output.cols()
                     << " frames.";
  }
  return true;
}

}  // namespace audition
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
