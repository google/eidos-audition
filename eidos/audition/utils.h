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

#ifndef EIDOS_AUDITION_UTILS_H_
#define EIDOS_AUDITION_UTILS_H_

#include <string>
#include <vector>

#include "eidos/port/logging.h"

namespace eidos {
namespace audition {

// Generates a vector of <n> uniformly distributed random real numbers in the
// supplied range [min_val, max_val].
std::vector<double> GetRandomUniformVector(int n, double min_val,
                                           double max_val);

// Generates a vector of <n> real numbers drawn from standard normal
// distribution characterized by the supplied mean and standard deviation.
std::vector<double> GetRandomNormalVector(int n, double mu, double sigma);

// Checks that the value is within range.
template <typename T>
void CheckWithinRange(T value, T range_min, T range_max,
                      const std::string &error_message_prefix) {
  if (value < range_min || value > range_max) {
    GOOGLE_LOG(FATAL) << error_message_prefix << " (value is " << value
                      << " should be within [" << range_min << ", " << range_max
                      << "] range";
  }
}

}  // namespace audition
}  // namespace eidos

#endif  // EIDOS_AUDITION_UTILS_H_

// Local Variables:
// mode: c++
// End:
