# Copyright 2019 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Simple library for reading and writing wave files.

cc_library(
    name = "wave",
    srcs = [
        "wave/file.cc",
        "wave/header.cc",
        "wave/header/data_header.cc",
        "wave/header/fmt_header.cc",
        "wave/header/riff_header.cc",
        "wave/header/wave_header.cc",
        "wave/header_list.cc",
    ],
    hdrs = [
        "wave/error.h",
        "wave/file.h",
        "wave/header.h",
        "wave/header/data_header.h",
        "wave/header/fmt_header.h",
        "wave/header/riff_header.h",
        "wave/header/wave_header.h",
        "wave/header_list.h",
    ],
    visibility = ["//visibility:public"],
)

# Local Variables:
# mode: python
# End:
