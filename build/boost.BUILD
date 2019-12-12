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

# Boost provides free peer-reviewed portable C++ source libraries.
#
# The libraries work well with the C++ Standard Library. Boost libraries are
# intended to be widely useful, and usable across a broad spectrum of
# applications. The Boost license encourages both commercial and
# non-commercial use.

load("@com_google_eidos_audition//build:boost.bzl", "boost_library")

cc_library(
    name = "headers",
    hdrs = glob([
        "boost/**/*.h",
        "boost/**/*.hpp",
        "boost/**/*.ipp",
    ]),
    # includes: list of include dirs added to the compile line
    includes = [
        ".",
    ],
    visibility = ["//visibility:public"],
)

boost_library(
    name = "filesystem",
    deps = [
        ":system",
    ],
)

boost_library(
    name = "graph",
    deps = [
        ":headers",
    ],
)

boost_library(
    name = "program_options",
    deps = [
        ":headers",
    ],
)

boost_library(
    name = "regex",
    deps = [
        ":headers",
    ],
)

boost_library(
    name = "system",
    deps = [
        ":headers",
    ],
)

boost_library(
    name = "thread",
    extra_hdrs = [
        "libs/thread/src/pthread/once_atomic.cpp",
    ],
    # Add source files for the pthread backend
    extra_srcs = glob([
        "libs/thread/src/pthread/once.cpp",
        "libs/thread/src/pthread/thread.cpp",
    ]),
    linkopts = [
        "-pthread",
    ],
    deps = [
        ":headers",
    ],
)

# Local Variables:
# mode: python
# End:
