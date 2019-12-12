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

#include "eidos/audition/psychoacoustic_scale_utils.h"

#include <cmath>

namespace eidos {
namespace audition {
namespace erb {

// This equation is an approximation for the moderate sound levels, normally
// between 10 Hz and 10 kHz. This equation is approximation of an earlier
// proposed square polynomial model ERB(f) = a * f^2 + b * f + c, originally
// proposed by Moore and Glasberg (1983).
double HzToErb(double hz) {
  return 24.7f * (0.00437 * hz + 1.0);
}

double ErbToHz(double erb) {
  return (erb / 24.7 - 1.0) / 0.00437;
}

double HzToErbRate(double hz) {
  return (21.4 * std::log10(hz * 0.00437 + 1.0));
}

double ErbRateToHz(double erb_rate) {
  return ((std::pow(10.0, erb_rate / 21.4) - 1.0) / 0.00437);
}

std::vector<double> GetFrequenciesHz(int num_channels,
                                     double min_frequency,
                                     double max_frequency) {
  const double min_erb_rate = erb::HzToErbRate(min_frequency);
  const double max_erb_rate = erb::HzToErbRate(max_frequency);
  const double step = (max_erb_rate - min_erb_rate) / num_channels;
  double erb_rate = min_erb_rate + step / 2.0;
  std::vector<double> center_frequencies(num_channels);
  for (int i = 0; i < num_channels; ++i) {
    center_frequencies[i] = erb::ErbRateToHz(erb_rate);
    erb_rate += step;
  }
  return center_frequencies;
}

}  // namespace erb
}  // namespace audition
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
