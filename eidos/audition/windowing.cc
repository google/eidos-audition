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

#include "eidos/audition/windowing.h"

#include <cmath>

#include "eidos/stubs/logging.h"

namespace eidos {
namespace audition {
namespace {

// Computes number of frames remaining in the waveform starting from the current
// offset.
int NumFramesInWave(int total_num_samples, int cur_offset, int frame_size,
                    int frame_shift) {
  if (frame_size + cur_offset > total_num_samples) {
    return 0;
  }
  return (total_num_samples - cur_offset - frame_size) / frame_shift + 1;
}

// Computes Hann window with dimension (num_channels, frame_size).
Eigen::ArrayXXd ComputeWindowHann(int frame_size, int num_channels) {
  const int n = frame_size - 1;
  Eigen::VectorXd vec(frame_size);
  for (int i = 0; i < frame_size; ++i) {
    vec(i) = 0.5 - 0.5 * std::cos(2.0 * M_PI * i / n);
  }
  return vec.replicate(/* row factor */1, num_channels).transpose();
}

}  // namespace

Eigen::ArrayXXd ComputeWindowFunction(WindowFunction window_function,
                                      int frame_size, int num_channels) {
  switch (window_function) {
    case WINDOW_FUNCTION_HANN:
      return ComputeWindowHann(frame_size, num_channels);
    case WINDOW_FUNCTION_NONE:
    default:
      return Eigen::ArrayXXd::Constant(num_channels, frame_size, 1.0);
  }
}

FrameInfo GetFrameInfo(const Eigen::ArrayXXd &input,
                       const StimulusConfig &config) {
  const int sample_rate = config.sample_rate();
  const double window_duration_sec = config.window_duration_sec();
  const double frame_shift_sec = config.frame_shift_sec();

  FrameInfo info;
  const double sample_period = 1.0 / sample_rate;  // In seconds.
  info.frame_size = static_cast<int>(window_duration_sec / sample_period);
  info.frame_shift = static_cast<int>(frame_shift_sec / sample_period);
  info.num_frames = NumFramesInWave(input.cols(), /* cur_offset */0,
                                    info.frame_size, info.frame_shift);
  return info;
}

std::vector<Eigen::ArrayXXd> Window(const Eigen::ArrayXXd &input,
                                    const StimulusConfig &config) {
  const FrameInfo &info = GetFrameInfo(input, config);
  GOOGLE_CHECK_GT(info.num_frames, 0)
      << "Invalid number of frames: " << info.num_frames;
  const int num_channels = input.rows();
  const Eigen::ArrayXXd &window_function = ComputeWindowFunction(
      config.window_function(), info.frame_size, num_channels);
  std::vector<Eigen::ArrayXXd> output;
  output.resize(info.num_frames);
  int cur_offset = 0, cur_frame_size = info.frame_size;
  const int num_samples = input.cols();
  for (int i = 0; i < info.num_frames; ++i) {
    if (cur_offset + info.frame_size > num_samples) {
      // Copy the remaining buffer and pad with zeros.
      cur_frame_size = num_samples - cur_offset;
      output[i] = Eigen::ArrayXXd::Zero(num_channels, info.frame_size);
      output[i].block(0, 0, num_channels, cur_frame_size) =
          input.block(0, cur_offset, num_channels, cur_frame_size);
    } else {
      output[i] = input.block(0, cur_offset, num_channels, cur_frame_size);
    }
    output[i] *= window_function;
    cur_offset += info.frame_shift;
  }
  return output;
}

Eigen::ArrayXXd WindowAndIntegrateTime(const Eigen::ArrayXXd &input,
                                       const StimulusConfig &config) {
  const std::vector<Eigen::ArrayXXd> &frames = Window(input, config);
  const int num_frames = frames.size();
  GOOGLE_CHECK_GT(num_frames, 0) << "Invalid number of frames: " << num_frames;
  const int frame_size = frames[0].cols();
  GOOGLE_CHECK_GT(frame_size, 0) << "Invalid frame size: " << frame_size;
  const int num_channels = frames[0].rows();
  Eigen::ArrayXXd output = Eigen::ArrayXXd(num_channels, num_frames);
  // TODO(agutkin): Row-wise norm (frames[i].rowwise().norm()) does not seem to
  // compile at the moment. Resorting to brute-force logic.
  for (int i = 0; i < num_frames; ++i) {
    const Eigen::ArrayXXd &frame = frames[i];
    for (int j = 0; j < num_channels; ++j) {
      double sum = 0.0;
      for (int k = 0; k < frame_size; ++k) {
        sum += (frame(j, k) * frame(j, k));
      }
      output(j, i) = std::sqrt(sum / frame_size);
    }
  }
  return output;
}

}  // namespace audition
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
