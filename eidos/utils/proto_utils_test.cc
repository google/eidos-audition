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

#include <unistd.h>
#include <filesystem>
#include <cstdio>
#include <string>
#include <utility>

#include "absl/strings/str_cat.h"
#include "eidos/stubs/status-matchers.h"
#include "eidos/utils/proto_utils.h"
#include "eidos/utils/test.pb.h"
#include "google/protobuf/util/message_differencer.h"
#include "gtest/gtest.h"

using google::protobuf::util::MessageDifferencer;

namespace eidos {
namespace utils {
namespace {

const char *const kTestDir = "eidos/utils/testdata";
const char *const kProtoName = "123";
const char *const kInvalidPath = "invalid";
constexpr int kProtoId = 456;

void SetupProto(TestProto *proto) {
  proto->set_name(kProtoName);
  proto->set_id(kProtoId);
}

}  // namespace

TEST(ProtoUtilsTest, ReadingTextFormat) {
  TestProto proto;
  SetupProto(&proto);
  const std::string proto_path = absl::StrCat(kTestDir, "/test.textproto");
  TestProto new_proto;
  EXPECT_FALSE(ReadTextProto(kInvalidPath, &new_proto).ok());
  ASSERT_OK(ReadTextProto(proto_path, &new_proto));
  EXPECT_EQ(new_proto.name(), proto.name());
  EXPECT_EQ(new_proto.id(), proto.id());
}

TEST(ProtoUtilsTest, WriteBinaryFormat) {
  TestProto proto;
  SetupProto(&proto);
  const std::string temp_path = std::tmpnam(nullptr);
  ASSERT_OK(WriteBinaryProto(proto, temp_path));
  TestProto new_proto;
  EXPECT_OK(ReadBinaryProto(temp_path, &new_proto));
  EXPECT_TRUE(MessageDifferencer::Equals(new_proto, proto));
  EXPECT_TRUE(std::filesystem::remove(temp_path));
}

TEST(ProtoUtilsTest, WriteTextFormat) {
  TestProto proto;
  SetupProto(&proto);
  const std::string temp_path = std::tmpnam(nullptr);
  ASSERT_OK(WriteTextProto(proto, temp_path));
  TestProto new_proto;
  EXPECT_OK(ReadTextProto(temp_path, &new_proto));
  EXPECT_TRUE(MessageDifferencer::Equals(new_proto, proto));
  EXPECT_TRUE(std::filesystem::remove(temp_path));
}

}  // namespace utils
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
