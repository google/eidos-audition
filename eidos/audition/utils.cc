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

#include "eidos/audition/utils.h"

#include <algorithm>
#include <random>

namespace eidos {
namespace audition {

std::vector<double> GetRandomUniformVector(int n, double min_val,
                                           double max_val) {
  std::random_device rnd_seed;  // Obtain a seed for the random number engine.
  std::mt19937 rnd_gen(rnd_seed());  // Standard Mersenne Twister engine.
  std::uniform_real_distribution<double> dist(min_val, max_val);
  std::vector<double> vec(n);
  std::generate(vec.begin(), vec.end(),
                [&dist, &rnd_gen]() { return dist(rnd_gen); });
  return vec;
}

std::vector<double> GetRandomNormalVector(int n, double mu, double sigma) {
  std::random_device rnd_seed;  // Obtain a seed for the random number engine.
  std::mt19937 rnd_gen(rnd_seed());  // Standard Mersenne Twister engine.
  std::normal_distribution<double> dist(mu, sigma);
  std::vector<double> vec(n);
  std::generate(vec.begin(), vec.end(),
                [&dist, &rnd_gen]() { return dist(rnd_gen); });
  return vec;
}

}  // namespace audition
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
