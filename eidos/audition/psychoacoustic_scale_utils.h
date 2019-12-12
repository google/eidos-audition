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

// Miscellaneous utility functions for various hearing critical bands inspired
// by psychoacoustic research.

#ifndef EIDOS_AUDITION_PSYCHOACOUSTICS_SCALE_UTILS_H_
#define EIDOS_AUDITION_PSYCHOACOUSTICS_SCALE_UTILS_H_

#include <vector>

namespace eidos {
namespace audition {

// Equivalent Rectangular Bandwidth (ERB):
// ---------------------------------------
// On ERB please consult:
//   - B.C.J. Moore and B.R. Glasberg (1996), "A revision of Zwicker's loudness
//     model", Acta Acustica united with Acustica, vol. 82, no. 2:335-345.
//   - B.C.J. Moore and B.R. Glasberg (1983), "Suggested formulae for
//     calculating auditory-filter bandwidths and excitation patterns",
//     J. Acoust. Soc. Am. 74: 750-753.
//   - B.R. Glasberg and B.C.J. Moore (1990), "Derivation of auditory filter
//     shapes from notched-noise data", Hearing Research, 47: 103-138.
//
// For a nice overview of ERB please see Lin, Shoufeng. "Logarithmic Frequency
// Scaling and Consistent Frequency Coverage for the Selection of Auditory
// Filterbank Center Frequencies.", arXiv preprint arXiv:1801.00075 (2017).
namespace erb {

// Converts characteristic frequency in Hz to ERB.
double HzToErb(double hz);

// Converts ERB to characteristic frequency in Hz.
double ErbToHz(double erb);

// Converts Hertz to Equivalent Rectangular Bandwidth (ERB) rate or scale, where
// this quantity ERBS(f) is defined as the number of ERBs below the supplied
// frequency in Hertz. See Glasberg and Moore (1996).
double HzToErbRate(double hz);

// Converts ERB rate back to Hertz.
double ErbRateToHz(double erb_rate);

// Creates an array of <num_channels> center frequencies (in Hz) equally
// spaced along the ERB scale in ascending order between <min_frequency> and
// <max_frequency>.
std::vector<double> GetFrequenciesHz(int num_channels,
                                     double min_frequency,
                                     double max_frequency);

}  // namespace erb
}  // namespace audition
}  // namespace eidos

#endif  // EIDOS_AUDITION_PSYCHOACOUSTICS_SCALE_UTILS_H_

// Local Variables:
// mode: c++
// End:
