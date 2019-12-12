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

# Defines the libsamplerate library that can be used inside Bazel.

def libsamplerate_library():
    native.exports_files(["build_libsamplerate.sh"])

    libname = "samplerate"
    tool = "build_lib%s.sh" % libname
    native.genrule(
        name = "build_lib%s" % libname,
        outs = [
             "lib%s.a" % libname,
             "%s.h" % libname,
        ],
        cmd = "$(location :%s) $(@D)" % tool,
        tools = [":%s" % tool],
    )

    native.cc_library(
        name = "libsamplerate",
        srcs = ["lib%s.a" % libname],
        hdrs = ["%s.h" % libname],
    )

# Local Variables:
# mode: python
# End:
