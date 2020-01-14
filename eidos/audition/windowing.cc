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

}  // namespace

FrameInfo GetFrameInfo(const Eigen::ArrayXXd &input,
                       int sample_rate,
                       double window_duration,
                       double frame_shift_sec) {
  FrameInfo info;
  const double sample_period = 1.0 / sample_rate;  // In seconds.
  info.frame_size = static_cast<int>(window_duration / sample_period);
  info.frame_shift = static_cast<int>(frame_shift_sec / sample_period);
  info.num_frames = NumFramesInWave(input.cols(), /* cur_offset */0,
                                    info.frame_size, info.frame_shift);
  return info;
}

std::vector<Eigen::ArrayXXd> Window(const Eigen::ArrayXXd &input,
                                    int sample_rate,
                                    double window_duration,
                                    double frame_shift_sec) {
  const FrameInfo &info = GetFrameInfo(input, sample_rate, window_duration,
                                       frame_shift_sec);
  GOOGLE_CHECK_GT(info.num_frames, 0)
      << "Invalid number of frames: " << info.num_frames;
  std::vector<Eigen::ArrayXXd> output;
  output.resize(info.num_frames);
  int cur_offset = 0, cur_frame_size = info.frame_size;
  const int num_channels = input.rows();
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
    cur_offset += info.frame_shift;
  }
  return output;
}

}  // namespace audition
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
