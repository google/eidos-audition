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

// Various testing helpers.

#ifndef EIDOS_UTILS_TEST_UTILS_H_
#define EIDOS_UTILS_TEST_UTILS_H_

#include <string>
#include <utility>

namespace eidos {
namespace utils {

// Returns a temporary path as a pair, where the first element is a file
// descriptor and a second element is a path name.
std::pair<int, std::string> TempPath();

}  // namespace utils
}  // namespace eidos

#endif  // EIDOS_UTILS_TEST_UTILS_H_

// Local Variables:
// mode: c++
// End:
