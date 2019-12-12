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

"""End-to-end test for cochlear feature extractor."""

import numpy as np
import os
import subprocess
import tempfile
import unittest

# Path to the test waveform in the original sample rate (24 kHz).
_TEST_WAVEFORM = "eidos/audition/testdata/ne_np/nep_letter_ma.wav"

# Path to the upsampled test waveform (100 kHz).
_TEST_WAVEFORM_100kHz = (
  "eidos/audition/testdata/ne_np/nep_letter_ma_100kHz.wav")

# Path to the auditory feature extractor configuration.
_TEST_CONFIG_PROTO_PATH = (
  "eidos/audition/testdata/configs/feature_extractor_config.textproto")

# Path to the binary.
_TOOL_PATH = "eidos/audition/auditory_feature_extractor"

# Simple configuration proto as (configurable) string.
_TEST_CONFIG_PROTO_CONTENTS = (
  "config { num_channels: %d } pipeline { models { model_type: MODEL_%s } }")

class AuditoryFeatureExtractorE2eTest(unittest.TestCase):
  def setUp(self):
    self.assertTrue(os.path.isfile(_TEST_WAVEFORM), "Waveform not found")
    self.assertTrue(os.path.isfile(_TOOL_PATH), "Tool binary not found")

  def _output_npy(self, filterbank, filterbank_options,
                  num_expected_channels, high_sample_rate=True,
                  config_proto_file=False, config_proto_contents=False):
    """Generates auditory signal using the supplied filterbank.

    Args:
      filterbank: Type of the filterbank.
      filterbank_options: Extra options for the filterbank.
      num_expected_channels: Expected number of filterbank channels.
      high_sample_rate: Whether to use upsampled audio. Some filterbanks,
        such as CARFAC, can use native sample rate. The rest require high
        sampling rates.
      config_proto_file: If enabled, checks initialization from a proto file.
    """
    with tempfile.NamedTemporaryFile("w") as out_file:
      features_file = out_file.name
      cmd = [
        _TOOL_PATH,
        "--output_file_format", "npy",
        "--waveform_file", (_TEST_WAVEFORM_100kHz if high_sample_rate
                            else _TEST_WAVEFORM),
        "--output_file", features_file,
      ]
      if not config_proto_file and not config_proto_contents:
        cmd.extend(["--models", filterbank])
      elif config_proto_file:
        cmd.extend(["--config_proto_file", _TEST_CONFIG_PROTO_PATH])
      else:  # Specify the config proto as string.
        cmd.extend(["--config_proto_contents",
                    _TEST_CONFIG_PROTO_CONTENTS % (
                      num_expected_channels, filterbank)])

      cmd.extend(filterbank_options)
      result = subprocess.call(cmd)
      self.assertEquals(0, result)
      self.assertTrue(os.path.isfile(features_file))
      features = np.load(features_file)
      dims = np.shape(features)
      # Number of columns corresponds to the number of samples.
      self.assertLess(1000, dims[1])
      # Number of rows corresponds to the auditory filterbank features
      # extracted for the requested sample rate.
      self.assertEquals(num_expected_channels, dims[0])

  def test_output_npy_carfac(self):
    """Checks CARFAC filters."""
    self._output_npy("CARFAC",
                     [],
                     72,  # Filterbank channels for 24 kHz waveform.
                     high_sample_rate=False)  # Can use native sample rate

  def test_output_npy_baumgarte(self):
    """Checks basic functioning of Frank Baumgarte's model."""
    self._output_npy("BAUMGARTE", ["--downsample_step", "1"],
                     251)  # Number of channels.

  def test_output_npy_zilany_2014(self):
    """Checks Zilany, et. al. (2014) model."""
    self._output_npy("ZILANY_IHC_2014", [],
                     251)  # Number of channels.

  def test_config_proto_file(self):
    """Checks simple configuration specified via proto file in text format."""
    self._output_npy("",   # Pipeline is defined in the proto.
                     [],
                     # Default number of channels for the Gammatone filter
                     # used by the pipeline defined in the configuration proto.
                     251,
                     high_sample_rate=False,  # Can use native sample rate
                     config_proto_file=True)

  def test_config_proto_contents(self):
    """Checks simple configuration specified as text string."""
    self._output_npy("GAMMATONE_SLANEY", [],
                     100,  # Number of channels.
                     high_sample_rate=True,  # Can use native sample rate
                     config_proto_file=False, config_proto_contents=True)


if __name__ == "__main__":
    unittest.main()

# Local Variables:
# mode: Python
# End:
