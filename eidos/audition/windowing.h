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

#ifndef EIDOS_AUDITION_WINDOWING_H_
#define EIDOS_AUDITION_WINDOWING_H_

#include <vector>

#include "Eigen/Core"
#include "eidos/audition/stimulus_config.pb.h"

namespace eidos {
namespace audition {

// Signal frame information.
struct FrameInfo {
  // Frame size (in samples).
  int frame_size = -1;

  // Frame shift (in samples): Advance frame by specified number of samples.
  int frame_shift = -1;

  // Number of frames in a signal.
  int num_frames = -1;
};

// Computes the window function given the window function type, the window
// width and the number of channels. Returns a window with dimension
// (num_channels, frame_size).
Eigen::ArrayXXd ComputeWindowFunction(WindowFunction window_function,
                                      int frame_size, int num_channels);

// Computes the frame information given the input signal, the duration of the
// window and the frame shift, both given in seconds.
FrameInfo GetFrameInfo(const Eigen::ArrayXXd &input,
                       const StimulusConfig &config);

// Runs windowing function over a two-dimensional input signal with dimensions
// (num_channels, num_samples). Returns a vector of <num_frames> two-dimensional
// arrays with dimension (num_channels, frame_size). The window width (in
// seconds) is provided by <window_duration_sec>. The window is advanced in
// steps specified by <frame_shift_sec>, specified in seconds.
std::vector<Eigen::ArrayXXd> Window(const Eigen::ArrayXXd &input,
                                    const StimulusConfig &config);

// Same as above, but also integrates across time dimension in each window by
// computing RMS for each channel. Returns an array with dimensions
// (num_channels, num_frames).
Eigen::ArrayXXd WindowAndIntegrateTime(const Eigen::ArrayXXd &input,
                                       const StimulusConfig &config);

}  // namespace audition
}  // namespace eidos

#endif  // EIDOS_AUDITION_WINDOWING_H_

// Local Variables:
// mode: c++
// End:
