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

// Utilities for manipulating the protos.

#ifndef EIDOS_UTILS_PROTO_UTILS_H_
#define EIDOS_UTILS_PROTO_UTILS_H_

#include <string>

#include "absl/strings/string_view.h"
#include "eidos/stubs/logging.h"
#include "google/protobuf/message.h"
#include "google/protobuf/stubs/status.h"
#include "google/protobuf/text_format.h"

namespace eidos {
namespace utils {

// Reads contents of a file into a string buffer. If <binary_mode> is enabled,
// opens the file in binary mode.
bool ReadFileContents(absl::string_view filepath, bool binary_mode,
                      std::string *contents);

// Writes contents to a file. If <binary_mode> is enabled the file is written
// in binary format.
bool WriteFileContents(const std::string &contents, bool binary_mode,
                       absl::string_view filepath);

// Reads binary protocol buffer message of the supplied type from a given file
// path.
template <typename ProtoType>
bool ReadBinaryProto(absl::string_view filepath, ProtoType *message) {
  // TODO(agutkin): It may be more efficient to parse directly from
  // std::istream.
  std::string contents;
  if (!ReadFileContents(filepath, true /* binary_mode */, &contents)) {
    return false;
  }
  return message->ParseFromString(contents);
}

// Reads textual protocol buffer message of the supplied type from a given file
// path.
template <typename ProtoType>
bool ReadTextProto(absl::string_view filepath, ProtoType *message) {
  // TODO(agutkin): It may be more efficient to parse directly from
  // io::ZeroCopyInputStream rather than from a string buffer.
  std::string contents;
  if (!ReadFileContents(filepath, false /* binary_mode */, &contents)) {
    return false;
  }
  using google::protobuf::TextFormat;
  return TextFormat::ParseFromString(contents, message);
}

// Writes proto to a given location in binary format.
template <typename ProtoType>
bool WriteBinaryProto(const ProtoType &message, absl::string_view filepath) {
  std::string contents;
  if (!message.SerializeToString(&contents)) {
    GOOGLE_LOG(ERROR) << "Binary serialization failed!";
    return false;
  }
  return WriteFileContents(contents, true /* binary mode */, filepath);
}

// Writes proto to a given location in text format.
template <typename ProtoType>
bool WriteTextProto(const ProtoType &message, absl::string_view filepath) {
  using google::protobuf::TextFormat;
  std::string contents;
  if (!TextFormat::PrintToString(message, &contents)) {
    GOOGLE_LOG(ERROR) << "Text serialization failed!";
    return false;
  }
  return WriteFileContents(contents, false /* binary mode */, filepath);
}

}  // namespace utils
}  // namespace eidos

#endif  // EIDOS_UTILS_PROTO_UTILS_H_

// Local Variables:
// mode: c++
// End:
