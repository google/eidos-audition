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

// Simple utility for extracting cochlear features.
//
// Example usage:
// --------------
//   - Using the model names:
//       bazel-bin/eidos/audition/auditory_feature_extractor \
//         --models GAMMATONE_SLANEY,SUMNER_SYNAPSE_2002 \
//         --waveform_file input.wav --output_file output.npy
//   - Using pipeline configuration contents:
//       bazel-bin/eidos/audition/auditory_feature_extractor \
//         --config_proto_contents "pipeline { models { model_type: MODEL_BAUMGARTE } }" \
//         --waveform_file input.wav --output_file output.npy
//   - Using pipeline configuration file:
//       bazel-bin/eidos/audition/auditory_feature_extractor \
//         --config_proto_file pipeline.textproto \
//         --waveform_file input.wav \
//         --output_file_format npz --output_file output.npz
//   - Processing multiple audio files:
//       bazel-bin/eidos/audition/auditory_feature_extractor \
//         --models BAUMGARTE,SUMNER_SYNAPSE_2002 \
//         --waveform_file_list file_list.txt --output_dir /tmp \
//         --output_file_format npz

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "absl/strings/ascii.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/str_split.h"
#include "eidos/audition/auditory_feature_extractor.h"
#include "eidos/audition/auditory_model_config.pb.h"
#include "eidos/audition/auditory_pipeline_config.pb.h"
#include "eidos/audition/auditory_response.h"
#include "eidos/audition/waveform.pb.h"
#include "eidos/audition/waveform_utils.h"
#include "eidos/stubs/logging.h"
#include "eidos/utils/proto_utils.h"
#include "google/protobuf/text_format.h"

ABSL_FLAG(std::string, waveform_file, "",
          "Audio file in WAVE format which is a subset of RIFF specification. "
          "Multiple sampling rates are supported but some models may require "
          "higher sampling rates, so the input must be upsampled. Number of "
          "bits per sample is not restricted to 16 bits. Only single channel "
          "audio is supported at the moment.");

ABSL_FLAG(std::string, waveform_file_list, "",
          "A text file containing newline-separated list of full paths to "
          "individual audio files. Mutually exclusive with --waveform_file "
          "flag.");

ABSL_FLAG(int, downsample_step, -1,
          "Simple temporal downsampling: Output every N-th sample. Please "
          "note, this flag is mutually exclusive with the proper resampling "
          "flags --output_resample_up_factor and "
          "--output_resample_down_factor.");

ABSL_FLAG(std::string, models, "",
          "Comma-separated list of model types. Supported types:\n"
          "\"MOCK_BASILAR_MEMBRANE\", \"CARFAC\", \"BAUMGARTE\", "
          "\"ZILANY_IHC_2014\", \"GAMMATONE_SLANEY\", "
          "\"MEDDIS_SYNAPSE_1986\", \"BRUCE_SYNAPSE_2018\", "
          "\"SUMNER_SYNAPSE_2002\", \"BRUCE_SPIKES_2018\", "
          "\"JACKSON_SPIKES\", \"ZHANG_SPIKES_2001\".");

ABSL_FLAG(std::string, config_proto_file, "",
          "Feature extractor configuration file in textual proto format.");

ABSL_FLAG(std::string, config_proto_contents, "",
          "Feature extractor configuration string as textual proto.");

ABSL_FLAG(double, audio_scaling_gain, 0.07,
          "Gain factor for audio scaling. Range: [0.01,0.1].");

ABSL_FLAG(std::string, output_file_format, "npy",
          "Format of the output file containing the model's response to "
          "stimulus. Supported values:\n"
          "\"npy\": Numpy .npy format: Single matrix with single-precision "
          "floats),\n"
          "\"npz\": Numpy .npz format: Multiple matrices with single-precision "
          "floats. Each matrix corresponds to the output of a particular "
          "auditory stage.");

ABSL_FLAG(std::string, output_file, "",
          "File path containing the output from the model in response to the "
          "individual stimulus. The format of the file is specified using "
          "--output_file_format option.");

ABSL_FLAG(std::string, output_dir, "",
          "Path to the directory where the model responses to individual "
          "stimulii provided by --waveform_file_list flag will be stored. "
          "Mutually exclusive with --output_file flag.");

ABSL_FLAG(int, num_channels, 251,
          "Number of channels (frequency bands) for analysis corresponding "
          "to <n> equidistant locations along the cochlea.");

ABSL_FLAG(double, min_cf_hz, 0.0,
          "Lowest characteristic frequency (CF) for analysis (in Hz). If "
          "unset, model-specific defaults will be used.");

ABSL_FLAG(double, max_cf_hz, 0.0,
          "Highest characteristic frequency (CF) for analysis (in Hz). If "
          "unset, model-specific defaults will be used.");

ABSL_FLAG(std::string, stimulus_config_proto_file, "",
          ("Text proto file specifying stimulus configuration. Can be used "
           "as an alternative to command-line flags."));

ABSL_FLAG(int, num_tasks, -1,
          "Number of tasks to run in parallel when processing a list of audio "
          "files. The value of -1 means that parallel processing is "
          "disabled. Set this flag to a non-negative number to enable the "
          "parallelism. A value of 0 will use the default hardware concurrency "
          "settings, positive number will force the precise number of tasks. "
          "Please note: Not all the models are thread-safe yet!");

ABSL_FLAG(bool, apply_window_to_outputs, false,
          "If enabled, applies windowing function to the response.");

ABSL_FLAG(double, window_duration_sec, 25E-3,
          "Window (also frame) duration (width) in seconds. Default is 25 ms.");

ABSL_FLAG(double, frame_shift_sec, 1E-2,
          "Frame shift (in seconds). After computing each frame, advance "
          "to the next by this amount. Default is 10 ms.");

ABSL_FLAG(std::string, window_function, "NONE",
          "Window function type specified as string. By default, no windowing "
          "function is applied. Supported values: \"NONE\", \"HANN\", "
          "\"HAMMING\".");

ABSL_FLAG(int, output_resample_up_factor, 0,
          "Resample the response: Upsampling factor <p>. The signal is "
          "resampled by <p/q>, where <q> is the downsampling factor.");

ABSL_FLAG(int, output_resample_down_factor, 0,
          "Resample the response: Downsampling factor <q>. The signal is "
          "resampled by <p/q>, where <p> is the upsampling factor.");

namespace eidos {
namespace audition {
namespace {

// Converts model type option string to the actual model types.
std::vector<AuditoryModelType> GetModelTypes() {
  const std::vector<std::string> model_strings = absl::StrSplit(
      absl::GetFlag(FLAGS_models), ',', absl::SkipEmpty());
  GOOGLE_CHECK(!model_strings.empty()) << "No models were supplied!";
  std::vector<AuditoryModelType> model_types;
  model_types.reserve(model_strings.size());
  AuditoryModelType model_type;
  for (const auto &model_string : model_strings) {
    GOOGLE_CHECK(AuditoryModelType_Parse(absl::StrCat(
        "MODEL_", model_string), &model_type))
        << "Failed to parse model type: " << model_string;
    model_types.emplace_back(model_type);
  }
  return model_types;
}

// Gets window function type from string.
WindowFunction GetWindowFunctionFromString() {
  const std::string function_name = absl::GetFlag(FLAGS_window_function);
  WindowFunction window = WINDOW_FUNCTION_NONE;
  GOOGLE_CHECK(WindowFunction_Parse(absl::StrCat(
      "WINDOW_FUNCTION_", function_name), &window))
      << "Failed to parse window function: " << function_name;
  return window;
}

// Initializes stimulus configuration from the command-line flags.
void InitStimulusConfig(StimulusConfig *stimulus_config) {
  const std::string stimulus_config_proto_file = absl::GetFlag(
      FLAGS_stimulus_config_proto_file);
  if (stimulus_config_proto_file.empty()) {
    stimulus_config->set_num_channels(absl::GetFlag(FLAGS_num_channels));
    stimulus_config->set_lowest_cf_hz(absl::GetFlag(FLAGS_min_cf_hz));
    stimulus_config->set_highest_cf_hz(absl::GetFlag(FLAGS_max_cf_hz));
    stimulus_config->set_audio_scaling_gain(absl::GetFlag(
        FLAGS_audio_scaling_gain));
    stimulus_config->set_downsample_step(absl::GetFlag(FLAGS_downsample_step));
    stimulus_config->set_apply_window_to_outputs(
        absl::GetFlag(FLAGS_apply_window_to_outputs));
    stimulus_config->set_window_duration_sec(
        absl::GetFlag(FLAGS_window_duration_sec));
    stimulus_config->set_frame_shift_sec(absl::GetFlag(FLAGS_frame_shift_sec));
    stimulus_config->set_window_function(GetWindowFunctionFromString());
    stimulus_config->set_output_resample_up_factor(
        absl::GetFlag(FLAGS_output_resample_up_factor));
    stimulus_config->set_output_resample_down_factor(
        absl::GetFlag(FLAGS_output_resample_down_factor));
  } else {
    GOOGLE_LOG(INFO) << "Reading stimulus configuration from "
                     << stimulus_config_proto_file;
    GOOGLE_CHECK(utils::ReadTextProto(stimulus_config_proto_file,
                                      stimulus_config))
        << "Failed to read stimulus configuration from: "
        << stimulus_config_proto_file;
  }
}

// Returns a fully qualified output path for the waveform input path given the
// output directory and a file format.
std::string OutputPathForInputWaveform(const std::string &input_wave_path,
                                       const std::string &output_file_format,
                                       const std::string &output_dir) {
  const std::filesystem::path path(input_wave_path);
  GOOGLE_CHECK(path.is_absolute()) << "Path should be absolute: "
                                   << input_wave_path;
  std::string output_path = absl::StrCat(output_dir, "/",
                                         path.stem().string());
  if (path.has_extension()) {
    output_path = absl::StrCat(output_path, ".", output_file_format);
  }
  return output_path;
}

}  // namespace

// Processes the input audio stimulus storing the result in the output file. The
// format of the output file is given by "output_file_format".
void ProcessFile(const std::string &waveform_file,
                 const std::string &output_file_format,
                 const std::string &output_file) {
  // Read the stimuli, denormalize the audio.
  GOOGLE_LOG(INFO) << "Reading " << waveform_file << " ...";
  Waveform waveform;
  GOOGLE_CHECK(ReadWaveform(waveform_file, &waveform))
      << "Failed to read waveform from \"" << waveform_file << "\"";
  DenormalizeWaveform(&waveform);

  // Fetch some of the options.
  const int num_channels = absl::GetFlag(FLAGS_num_channels);
  const std::string config_proto_file = absl::GetFlag(FLAGS_config_proto_file);
  const std::string config_proto_contents = absl::GetFlag(
      FLAGS_config_proto_contents);
  const bool npz_output_format = (output_file_format == "npz");

  // Prepare configuration.
  AuditoryFeatureExtractorConfig config;
  GOOGLE_CHECK_LT(0, num_channels)
      << "Invalid number of channels: " << num_channels;
  StimulusConfig *stimulus_config = config.mutable_config();
  InitStimulusConfig(stimulus_config);
  if (!absl::GetFlag(FLAGS_models).empty()) {
    // Configuration is specified via list of models.
    const std::vector<AuditoryModelType> &model_types = GetModelTypes();
    bool store_intermediate_outputs = false;
    if (model_types.size() > 1 && npz_output_format) {
      store_intermediate_outputs = true;
    }
    AuditoryPipelineConfig *pipeline_config = config.mutable_pipeline();
    pipeline_config->mutable_models()->Reserve(model_types.size());
    for (const auto model_type : model_types) {
      AuditoryModelConfig *model_config = pipeline_config->add_models();
      model_config->set_model_type(model_type);
    }
    stimulus_config->set_store_intermediate_outputs(store_intermediate_outputs);
  } else if (!config_proto_file.empty() || !config_proto_contents.empty()) {
    // Please note, the proto may override some of the stimulus configuration
    // options provided through the flags. Therefore we simply merge the
    // contents of the default config with the proto loaded through the pipeline
    // configuration flag.
    AuditoryFeatureExtractorConfig pipeline_config;
    if (!config_proto_file.empty()) {
      GOOGLE_LOG(INFO) << "Reading configuration from " << config_proto_file;
      GOOGLE_CHECK(utils::ReadTextProto(config_proto_file, &pipeline_config))
          << "Failed to read configuration from: " << config_proto_file;
    } else {  // Read from string.
      using google::protobuf::TextFormat;
      GOOGLE_CHECK(TextFormat::ParseFromString(config_proto_contents,
                                               &pipeline_config))
          << "Failed to parse config contents";
    }
    config.MergeFrom(pipeline_config);
    if (config.pipeline().models().size() > 1 && npz_output_format) {
      stimulus_config->set_store_intermediate_outputs(true);
    }
  } else {
    GOOGLE_LOG(FATAL) << "Specify --models, --config_proto_file or "
                      << "--config_proto_contents";
  }

  // Generate cochlear features.
  GOOGLE_LOG(INFO) << "Extracting features ...";
  AuditoryFeatureExtractor extractor;
  AuditoryResponse response;
  GOOGLE_CHECK(extractor.ComputeFeatures(waveform, config, &response))
      << "Failed to extract features from " << waveform_file;
  GOOGLE_CHECK_LE(1, response.outputs().size())
      << "Expected at least one output from the auditory pipeline";

  // Write the results.
  const std::string output_format_msg =
      (npz_output_format ? "(npz)" : "(npy)");
  GOOGLE_LOG(INFO) << "Saving " << output_format_msg << " output to "
                   << output_file << " ...";
  if (!npz_output_format) {
    response.WriteNpy(output_file);
  } else {  // Assume .npz output format.
    response.WriteNpz(output_file);
  }
}

void Run() {
  // Figure out the output file format.
  const std::string output_file_format = absl::GetFlag(
      FLAGS_output_file_format);
  GOOGLE_CHECK(output_file_format == "npy" || output_file_format == "npz")
      << "Unsupported output file format: " << output_file_format;

  // Process individual file or a list of files.
  const std::string waveform_file = absl::GetFlag(FLAGS_waveform_file);
  const std::string output_file = absl::GetFlag(FLAGS_output_file);
  const std::string waveform_file_list = absl::GetFlag(
      FLAGS_waveform_file_list);
  const std::string output_dir = absl::GetFlag(FLAGS_output_dir);
  if (!waveform_file.empty() && !output_file.empty()) {
    ProcessFile(waveform_file, output_file_format, output_file);
  } else if (!waveform_file_list.empty() && !output_dir.empty()) {
    // Check that output directory is valid.
    const std::filesystem::path output_dir_path(output_dir);
    GOOGLE_CHECK(std::filesystem::is_directory(output_dir_path))
        << "Output directory does not exist or is not a directory: "
        << output_dir;

    // Read the file list.
    std::string contents;
    GOOGLE_CHECK(utils::ReadFileContents(
        waveform_file_list, /* binary_mode= */false, &contents))
        << "Failed to read " << waveform_file_list;
    std::vector<std::string> waveform_files = absl::StrSplit(
        contents, '\n', absl::SkipEmpty());
    GOOGLE_CHECK(!waveform_files.empty()) << "No waveforms were supplied!";
    for (auto &wav_file : waveform_files) {
      absl::StripAsciiWhitespace(&wav_file);
    }
    std::sort(waveform_files.begin(), waveform_files.end());
    std::unique(waveform_files.begin(), waveform_files.end());

    // Process individual files.
    int num_parallel_tasks = absl::GetFlag(FLAGS_num_tasks);
    if (num_parallel_tasks < 0) {  // Sequential mode.
      for (const auto &wav_file : waveform_files) {
        ProcessFile(wav_file, output_file_format, OutputPathForInputWaveform(
            wav_file, output_file_format, output_dir));
      }
    } else {  // Parallel mode of operation.
      GOOGLE_LOG(WARNING) << "NOTE: Not all the models are thread-safe yet. "
                          << "Please check the code to see whether parallel "
                          << "mode is suitable for a particular model.";
      const int hardware_concurrency = std::thread::hardware_concurrency();
      if (num_parallel_tasks == 0 ||
          num_parallel_tasks > hardware_concurrency) {
        num_parallel_tasks = hardware_concurrency;
      }
      const int num_waveforms = waveform_files.size();
      int num_batches = num_waveforms / num_parallel_tasks;
      if (num_batches == 0) num_batches = 1;
      GOOGLE_LOG(INFO) << "Number of tasks: " << num_parallel_tasks
                       << ", batches: " << num_batches;

      // Process the waveforms in parallel using the minimal threading interface
      // of C++11. Ideally a thread pool will fit here, but dragging in
      // additional dependencies was ruled out.
      for (int i = 0; i < num_waveforms; ) {
        const int num_remaining = num_waveforms - i;
        if (num_remaining == 0) {
          break;
        } else if (num_parallel_tasks > num_remaining) {
          num_parallel_tasks = num_remaining;
        }
        std::vector<std::thread> workers(num_parallel_tasks);
        for (int j = 0; j < num_parallel_tasks; ++j) {
          const std::string &wav_file = waveform_files[i + j];
          const std::string &output_path = OutputPathForInputWaveform(
              wav_file, output_file_format, output_dir);
          workers[j] = std::thread(ProcessFile, wav_file, output_file_format,
                                   output_path);
        }
        for (auto &worker : workers) {
          if (worker.joinable()) worker.join();
        }
        i += num_parallel_tasks;
      }
    }
  } else {
    GOOGLE_LOG(FATAL) << "Wrong combination of input/output flags!";
  }
}

}  // namespace audition
}  // namespace eidos

int main(int argc, char *argv[]) {
  absl::ParseCommandLine(argc, argv);
  GOOGLE_CHECK(!absl::GetFlag(FLAGS_waveform_file).empty() ||
               !absl::GetFlag(FLAGS_waveform_file_list).empty())
      << "--waveform_file or --waveform_file_list are required";
  GOOGLE_CHECK(!absl::GetFlag(FLAGS_output_file).empty() ||
               !absl::GetFlag(FLAGS_output_dir).empty())
      << "Specify --output_file or --output_dir";
  eidos::audition::Run();
  return 0;
}

// Local Variables:
// mode: c++
// End:
