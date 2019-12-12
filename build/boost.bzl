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

DEFAULT_COPTS = ["-Wno-unused-value"]

def boost_library(name, defines = None, deps = None, extra_srcs = None, extra_hdrs = None, linkopts = None):
    if defines == None:
        defines = []

    if deps == None:
        deps = []

    if extra_srcs == None:
        extra_srcs = []

    if extra_hdrs == None:
        extra_hdrs = []

    if linkopts == None:
        linkopts = []

    return native.cc_library(
        name = name,
        visibility = ["//visibility:public"],
        defines = defines,
        hdrs = native.glob([
            x % name
            for x in [
                "libs/%s/src/*.hpp",
                "boost/%s/**/*.ipp",
            ]
        ]) + extra_hdrs,
        srcs = native.glob([
            x % name
            for x in [
                "libs/%s/src/*.cpp",
                "boost/%s/detail/*.hpp",
            ]
        ]) + extra_srcs,
        deps = deps,
        copts = DEFAULT_COPTS,
        linkopts = linkopts,
    )

# Local Variables:
# mode: python
# End:
