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

#include "eidos/utils/proto_utils.h"

#include <fstream>
#include <sstream>

#include "absl/strings/str_cat.h"

namespace eidos {
namespace utils {

absl::Status ReadFileContents(absl::string_view filepath, bool binary_mode,
                              std::string *contents) {
  const std::string filepath_str(filepath);
  std::ios_base::openmode mode = std::ios_base::in;
  if (binary_mode) mode |= std::ios::binary;
  std::ifstream input_stream(filepath_str, mode);
  if (!input_stream) {
    return absl::NotFoundError(absl::StrCat(
        "Failed to open \"", filepath_str, "\""));
  }
  std::ostringstream string_buffer;
  string_buffer << input_stream.rdbuf();
  *contents = string_buffer.str();
  return absl::OkStatus();
}

bool WriteFileContents(const std::string &contents, bool binary_mode,
                       absl::string_view filepath) {
  const std::string filepath_str(filepath);
  std::ios_base::openmode mode = std::ios_base::out;
  if (binary_mode) mode |= std::ios::binary;
  std::ofstream output_stream(filepath_str, mode);
  if (!output_stream) {
    GOOGLE_LOG(ERROR) << "Failed to open \"" << filepath_str << "\"";
    return false;
  }
  output_stream.write(contents.c_str(), contents.size());
  if (!output_stream) {
    GOOGLE_LOG(ERROR) << "Failed to write to \"" << filepath_str << "\"";
    return false;
  }
  return true;
}

}  // namespace utils
}  // namespace eidos
