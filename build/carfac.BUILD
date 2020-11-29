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

# The CAR-FAC (cascade of asymmetric resonators with fast-acting compression).

cc_library(
    name = "carfac",
    srcs = [
        "cpp/binaural_sai.cc",
        "cpp/carfac.cc",
        "cpp/ear.cc",
        "cpp/sai.cc",
    ],
    hdrs = [
        "cpp/agc.h",
        "cpp/binaural_sai.h",
        "cpp/car.h",
        "cpp/carfac.h",
        "cpp/carfac_util.h",
        "cpp/common.h",
        "cpp/ear.h",
        "cpp/ihc.h",
        "cpp/sai.h",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "@com_gitlab_eigen//:eigen",
    ],
)

cc_library(
    name = "test_util",
    hdrs = ["cpp/test_util.h"],
)

cc_test(
    name = "carfac_test",
    size = "small",
    srcs = ["cpp/carfac_test.cc"],
    data = [
        "test_data/agc_test-audio.txt",
        "test_data/agc_test-matlab-bm1.txt",
        "test_data/agc_test-matlab-bm2.txt",
        "test_data/agc_test-matlab-nap1.txt",
        "test_data/agc_test-matlab-nap2.txt",
        "test_data/binaural_test-audio.txt",
        "test_data/binaural_test-matlab-bm1.txt",
        "test_data/binaural_test-matlab-bm2.txt",
        "test_data/binaural_test-matlab-nap1.txt",
        "test_data/binaural_test-matlab-nap2.txt",
        "test_data/ihc_just_hwr_test-audio.txt",
        "test_data/ihc_just_hwr_test-matlab-bm1.txt",
        "test_data/ihc_just_hwr_test-matlab-bm2.txt",
        "test_data/ihc_just_hwr_test-matlab-nap1.txt",
        "test_data/ihc_just_hwr_test-matlab-nap2.txt",
        "test_data/long_test-audio.txt",
        "test_data/long_test-matlab-bm1.txt",
        "test_data/long_test-matlab-bm2.txt",
        "test_data/long_test-matlab-nap1.txt",
        "test_data/long_test-matlab-nap2.txt",
    ],
    deps = [
        ":carfac",
        ":test_util",
        "@com_google_googletest//:gtest_main",
    ],
)

cc_test(
    name = "sai_test",
    size = "small",
    srcs = ["cpp/sai_test.cc"],
    data = [
        "test_data/binaural_test-matlab-nap1.txt",
        "test_data/binaural_test-matlab-sai1.txt",
    ],
    deps = [
        ":carfac",
        ":test_util",
        "@com_google_googletest//:gtest_main",
    ],
)

cc_test(
    name = "binaural_sai_test",
    size = "small",
    srcs = ["cpp/binaural_sai_test.cc"],
    deps = [
        ":carfac",
        ":test_util",
        "@com_google_googletest//:gtest_main",
    ],
)

# Local Variables:
# mode: python
# End:
