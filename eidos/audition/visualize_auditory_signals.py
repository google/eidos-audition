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

"""A simple visualizer for auditory periphery signals."""

import os
import matplotlib as mpl
if os.environ.get('DISPLAY', '') == '':
  print("Note: Comment out Agg display backend for this tool to work "
        "on MacOS X.")
  mpl.use('Agg')

import matplotlib.axes as axes
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np

from absl import app
from absl import flags
from absl import logging

flags.DEFINE_string(
    "input_signal_file", None,
    ("Input file containing a single signal matrix (.npy file) or an "
     "associative container (.npz file) of signal matrixes. Each signal is "
     "stored in a double matrix of dimension N by M, where N is the number of "
     "samples (or frames, if windowing is enabled) and M is the number of "
     "filterbank channels."))

flags.DEFINE_integer(
    "input_sample_rate", -1,
    "Sample rate of the input signal at which it was produced.")

flags.DEFINE_bool(
    "time_smooth", False,
    "By default the response is displayed in samples. If enabled, smooth the "
    "response the display the output in seconds.")

flags.DEFINE_bool(
    "color_log_norm", False,
    "Use log normalization for colors.")

flags.DEFINE_bool(
    "show_colorbar", True,
    "Display color bar.")

flags.DEFINE_integer(
    "x_axis_start", -1,
    "The region offset along the time axis.")

flags.DEFINE_string(
    "color_map", "jet",
    "Color map for displaying the image (e.g., \"bone\", \"gray_r\").")

flags.DEFINE_integer(
    "x_axis_end", -1,
    "The right limit of the region along the time axis.")

flags.DEFINE_string(
    "output_image_file", None,
    "If specified, will output the PNG image to a file.")

FLAGS = flags.FLAGS

# Duration of analysis frame for the time smoothing in seconds.
_TIME_SMOOTHING_FRAME_LENGTH_SEC = 0.005

# Name of the array holding the central frequencies.
_CENTRAL_FREQUENCIES_NAME = "CENTER_FREQUENCIES"


def _x_limits(signal):
  """Computes the x-axis region to be displayed.

  Args:
     signal: Input signal.

  Returns:
     A tuple consisting of left and right limits.
  """
  xlim_left = 0
  xlim_right = signal.shape[1] - 1
  if FLAGS.x_axis_end >= signal.shape[1]:
    logging.fatal("Invalid time region: [%d, %d]",
                  FLAGS.x_axis_start, FLAGS.x_axis_end)
  elif FLAGS.x_axis_end > 0 and FLAGS.x_axis_end <= FLAGS.x_axis_start:
    logging.fatal("Invalid right rimit: [%d, %d]",
                  FLAGS.x_axis_start, FLAGS.x_axis_end)
  if FLAGS.x_axis_start > 0 and FLAGS.x_axis_end > 0:
    xlim_left = FLAGS.x_axis_start
    xlim_right = FLAGS.x_axis_end
  elif FLAGS.x_axis_start > 0:
    xlim_left = FLAGS.x_axis_start
  elif FLAGS.x_axis_end > 0:
    xlim_right = FLAGS.x_axis_end
  return xlim_left, xlim_right


def _time_smooth(signal):
  """Smoothes the input signal usually given in samples in time.

  In order to achieve this the input sample rate is required.

  Args:
     signal: Input signal in samples.

  Returns:
     Time-smoothed signal.
  """
  if FLAGS.input_sample_rate <= 0:
    logging.fatal("Invalid sampling rate: %d", FLAGS.input_sample_rate)
  samples_in_frame = int(
    _TIME_SMOOTHING_FRAME_LENGTH_SEC * FLAGS.input_sample_rate + 0.5)
  logging.info("Frame length: %d samples", samples_in_frame)
  num_frames = int(signal.shape[1] / samples_in_frame)
  logging.info("Duration: %d frames", num_frames)
  smoothed_signal = np.zeros(shape=(signal.shape[0], num_frames), dtype=np.float32)
  signal_offset = 0
  for n in range(num_frames):
    smoothed_signal[:, n] = np.average(
      signal[:, signal_offset:signal_offset + samples_in_frame], axis=1)
    signal_offset += samples_in_frame
  return smoothed_signal


def _signal_type_to_label(signal_type):
  """Converts signal type key to a readable string label.

  Args:
     signal_type: String key from an .npz file.

  Returns:
     Readable string.
  """
  if signal_type == 'OUTPUT_BASILAR_MEMBRANE_DISPLACEMENT':
    return "Basilar Membrane Displacement"
  elif signal_type == 'OUTPUT_IHC_TRANSMEMBRANE_POTENTIAL':
    return "IHC Transmembrane Potential"
  elif signal_type == 'OUTPUT_SYNAPSE_FIRING_RATES':
    return "Auditory Nerve Firing Rates"
  elif signal_type == 'OUTPUT_SYNAPSE_SPIKE_TIMES':
    return "Spike Times"
  else:
    logging.fatal("Unknown signal type: %s", signal_type)


def _display_signal(fig, subplot_index, num_subplots, signal, signal_type,
                    center_frequencies):
  """Displays the supplied signal.

  Args:
     fig: Figure where the plotting will happen.
     subplot_index: An integer specifying 1-based ID of the subplot.
     num_subplots: Number of subplots.
     signal: Two-dimensional matrix. The dimensions are number of channels by
       number of samples or frames.
     signal_type: String specifying the type of the signal. May be none if the
       data is loaded from an .npy file.
     center_frequencies: Vector containing center frequencies (CFs). May be
       empty if this information is missing from the input files.
  """
  # Prepare color normalization, if specified.
  color_norm = None
  if FLAGS.color_log_norm:
    vmax = abs(signal).max()
    vmin = vmax / 1.0E5
    color_norm = colors.LogNorm(vmin=vmin, vmax=vmax)

  # Prepare the figure.
  cmap = FLAGS.color_map
  origin = "lower"  # Vertical axis points upwards.
  ax = fig.add_subplot(num_subplots, 1, subplot_index)
  xlimits = _x_limits(signal)
  ax.set_xlim(left=xlimits[0], right=xlimits[1])
  pcm = ax.imshow(np.abs(signal), norm=color_norm, cmap=cmap, origin=origin,
                  aspect="auto")
  if signal_type:
    ax.set_title(_signal_type_to_label(signal_type))
  if subplot_index == num_subplots:
    # Only set the x-axis label for the last subplot.
    ax.set_xlabel("Time (%s)" % (
        "sec" if FLAGS.time_smooth else "samples or frames"))
  ax.set_ylabel("Channels")
  ax.autoscale(enable=False)
  pad = 0.15 if center_frequencies.size > 0 else 0.05
  fig.colorbar(pcm, ax=ax, pad=pad)
  if center_frequencies.size > 0:
    ax2 = ax.twinx()
    ax2.set_ylabel("Frequency (Hz)")
    ax2.set_ylim(bottom=center_frequencies[0], top=center_frequencies[-1])
    ax2.autoscale(enable=False)

  # Convert frames to seconds if smoothing was requested.
  if FLAGS.time_smooth:
    labels = ax.get_xticks().tolist()
    labels_sec = [str(
      int(item) * _TIME_SMOOTHING_FRAME_LENGTH_SEC) for item in labels]
    ax.set_xticklabels(labels_sec)


def main(unused_argv):
  # Read the signal from a numpy matrix and smooth it in time, if required.
  # The input signal is first converted to single-precision floats otherwise
  # the plotting library cannot handle double-precision and crashes.
  logging.info("Reading signal from %s ...", FLAGS.input_signal_file)
  signal_data = np.load(FLAGS.input_signal_file)
  signals = []
  center_frequencies = np.array([])
  if isinstance(signal_data, np.ndarray):
    signal = signal_data.astype('float32')
    logging.info("Read %s matrix.", signal.shape)
    if FLAGS.time_smooth:
      signal = _time_smooth(signal)
    signals.append((signal, None))
  else:  # Assume .npz.
    # Read center frequencies information, if present.
    if _CENTRAL_FREQUENCIES_NAME in signal_data.keys():
      center_frequencies = signal_data[_CENTRAL_FREQUENCIES_NAME]
      logging.info("Read %d center frequencies", center_frequencies.shape[0])

    # Process the signals.
    for signal_type in signal_data.keys():
      if signal_type == _CENTRAL_FREQUENCIES_NAME:
        continue
      signal = signal_data[signal_type]
      logging.info("%s: Read %s matrix.", signal_type, signal.shape)
      if FLAGS.time_smooth:
        signal = _time_smooth(signal)
      signals.append((signal, signal_type))

  # Plot everything.
  fig = plt.figure()
  logging.info("Displaying %d signals ...", len(signals))
  for subplot_i, (signal, signal_type) in enumerate(signals):
    logging.info("[%d] Plotting %s ...", subplot_i, signal_type)
    _display_signal(fig, subplot_i + 1, len(signals), signal, signal_type,
                    center_frequencies)

  # Save to a file or show the figure.
  plt.tight_layout()  # Makes sure the title doesn't overlap with x-axis ticks.
  if FLAGS.output_image_file:
    logging.info("Saving image to %s ...", FLAGS.output_image_file)
    fig.savefig(FLAGS.output_image_file)
  else:
    plt.show(block=True)


if __name__ == "__main__":
  flags.mark_flag_as_required("input_signal_file")
  app.run(main)


# Local Variables:
# mode: Python
# End:
