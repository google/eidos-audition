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

"""Sanity checks for auditory pipeline configurations."""

def config_tester(name, test_enabled = True):
    """Checks the supplied auditory feature extractor configuration.

    Args:
      name: Name for the test. This should match name of the existing pipeline
        configuration proto in the directory where the rule is loaded, excluding
        the extension.
      test_enabled: Flag indicating whether the test should be run.
    """
    if not test_enabled:
      print("%s: Configuration disabled." % name)
      return

    proto_filename = "%s.textproto" % name
    native.cc_test(
        name = "%s_smoke_test" % name,
        args = [
            "--config_proto_file",
            "$(location %s)" % proto_filename,
        ],
        data = [proto_filename],
        deps = ["//eidos/audition:auditory_feature_extractor_config_tester"],
    )


# Local Variables:
# mode: python
# End:
