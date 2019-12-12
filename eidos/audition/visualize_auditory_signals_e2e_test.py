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

"""Simple end-to-end test for signal visualizer."""

import matplotlib.image as mpimg
import os
import subprocess
import tempfile
import unittest

# Path to the test signal.
_TEST_SIGNAL = (
    "eidos/audition/testdata/auditory/1kHz_60dB_SPL_tone_baumgarte.npy")

# Path to the binary.
_TOOL_PATH = "eidos/audition/visualize_auditory_signals"


class VisualizeCochlearSignalsE2eTest(unittest.TestCase):
  def setUp(self):
    self.assertTrue(os.path.isfile(_TEST_SIGNAL), "Signal not found")
    self.assertTrue(os.path.isfile(_TOOL_PATH), "Tool binary not found")

  def test_image_generation(self):
    """Checks that we can convert the input signal to an image."""
    with tempfile.NamedTemporaryFile("w", suffix=".png") as out_file:
      image_file = out_file.name
      cmd = [
          _TOOL_PATH,
          "--input_signal_file", _TEST_SIGNAL,
          "--output_image_file", image_file,
          "--x_axis_start", str(50),
          "--x_axis_end", str(700),
      ]
      result = subprocess.call(cmd)
      self.assertEquals(0, result)
      self.assertTrue(os.path.isfile(image_file))
      image = mpimg.imread(image_file)
      self.assertLess(10, image.shape[0])
      self.assertLess(10, image.shape[1])


if __name__ == "__main__":
  unittest.main()

# Local Variables:
# mode: Python
# End:
