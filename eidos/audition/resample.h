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

#ifndef EIDOS_AUDITION_RESAMPLE_H_
#define EIDOS_AUDITION_RESAMPLE_H_

#include "Eigen/Core"
#include "eidos/audition/stimulus_config.pb.h"

namespace eidos {
namespace audition {

// Checks whether resampler should be applied based on the configuration.
// Please note, if <downsample_step> is bigger than one, this resampling
// API will not be called.
bool ShouldResampleOutputs(const StimulusConfig &config);

// Resamples the input signal with dimensions (num_channels, num_samples). The
// input sample rate and the resampling ratios  are provided by the
// configuration <config>. Returns the resampled signal with dimensions
// (num_cannels, new_num_samples), where each channel is resampled
// independently.
//
// The underlying resampler is from https://github.com/terrygta/SignalResampler.
// It reimplements Matlab one-dimensional <resample> function
//   http://www.mathworks.com.au/help/signal/ref/resample.html
// using open-source <upfirdn> implementation by Motorola
//   http://sourceforge.net/motorola/upfirdn/home/Home/.
//
// According to the original Matlab documentation:
//   y = resample(x,p,q) resamples the input sequence, x, at p/q times the
//   original sample rate. If x is a matrix, then resample treats each column
//   of x as an independent channel. resample applies an antialiasing FIR
//   lowpass filter to x and compensates for the delay introduced by the filter.
Eigen::ArrayXXd Resample(const Eigen::ArrayXXd &input,
                         const StimulusConfig &config);

}  // namespace audition
}  // namespace eidos

#endif  // EIDOS_AUDITION_RESAMPLE_H_

// Local Variables:
// mode: c++
// End:
