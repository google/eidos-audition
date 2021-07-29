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

#include "eidos/audition/resampler.h"

#include <vector>

#include "eidos/stubs/logging.h"
#include "resample.h"  // @com_github_resample

namespace eidos {
namespace audition {

bool ShouldResampleOutputs(const StimulusConfig &config) {
  return (
      (config.output_resample_up_factor() > 0) &&
      (config.output_resample_down_factor() > 0) &&
      (config.output_resample_up_factor() !=
       config.output_resample_down_factor()) &&
      (config.downsample_step() <= 1));
}

Eigen::ArrayXXd DownstepResample(const Eigen::ArrayXXd &input,
                                 const StimulusConfig &config) {
  const int downsample_step = config.downsample_step();
  if (downsample_step <= 1) {
    return input;
  }
  const int num_channels = input.rows();
  const int num_output_samples = input.cols() / downsample_step;
  Eigen::ArrayXXd output = Eigen::ArrayXXd(num_channels, num_output_samples);
  for (int i = 0; i < num_channels; ++i) {
    int n = 0;
    for (int j = 0; j < num_output_samples; ++j) {
      output(i, j) = input(i, n);
      n += downsample_step;
    }
  }
  return output;
}

Eigen::ArrayXXd Resample(const Eigen::ArrayXXd &input,
                         const StimulusConfig &config) {
  const int sample_rate = config.sample_rate();
  GOOGLE_CHECK_LT(0.0, sample_rate);
  const int up_factor = config.output_resample_up_factor();
  GOOGLE_CHECK_LT(0.0, up_factor) << "Invalid upsampling factor: " << up_factor;
  const int down_factor = config.output_resample_down_factor();
  GOOGLE_CHECK_LT(0.0, down_factor)
      << "Invalid downsampling factor: " << down_factor;
  if (up_factor == down_factor) {  // No-op.
    return input;
  }

  // Figure out the number of output samples by resampling the first channel.
  const int num_samples = input.cols();
  auto channel = std::vector<double>(input.row(0).data(),
                                     input.row(0).data() + num_samples);
  std::vector<double> output_channel;
  ::resample(up_factor, down_factor, channel, output_channel);
  const int num_output_samples = output_channel.size();
  GOOGLE_CHECK_LT(0, num_output_samples)
      << "Invalid number of output samples: " << num_output_samples;

  // Resample the signal starting from the second channel.
  const int num_channels = input.rows();
  Eigen::ArrayXXd output = Eigen::ArrayXXd(num_channels, num_output_samples);
  output.row(0) = Eigen::Map<Eigen::ArrayXd>(output_channel.data(),
                                             num_output_samples);
  for (int i = 1; i < num_channels; ++i) {
    auto channel = std::vector<double>(input.row(i).data(),
                                       input.row(i).data() + num_samples);
    ::resample(up_factor, down_factor, channel, output_channel);
    GOOGLE_CHECK_EQ(output_channel.size(), num_output_samples);
    output.row(i) = Eigen::Map<Eigen::ArrayXd>(output_channel.data(),
                                               num_output_samples);
  }
  return output;
}

}  // namespace audition
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
