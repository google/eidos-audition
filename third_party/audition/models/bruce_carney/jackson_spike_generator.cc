// Copyright 2003 B. Scott Jackson
// Copyright 2011 Muhammad S. A. Zilany
// Copyright 2019 Google LLC
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// version 2 as published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.

#include "third_party/audition/models/bruce_carney/jackson_spike_generator.h"

#include <algorithm>
#include <cstdint>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <random>

#include "eidos/audition/auditory_model_config.pb.h"
#include "eidos/port/logging.h"

namespace eidos {
namespace audition {
namespace {

// Dead time or absolute refractory period (sec).
constexpr double kDeadTime = 0.75E-3;

// Relative refractory parameters:
//
// According to Zhang and Carney (2001), following spike history parameters for
// two exponents fit
//   - Westerman and Smith (1989): "A diffusion model of the transient response
//     of the cochlear inner hair cell synapse."., JASA, vol. 83, no. 6.
//   - Gaumond and Kim (1983): "Response of cochlear nerve fibers to brief
//     acoustic stimuli: Role of discharge‐history effects.", JASA,
//     vol. 74, no. 5.
constexpr double kC0 = 0.5;     // Dimensionless.
constexpr double kS0 = 0.001;   // sec
constexpr double kC1 = 0.5;     // Dimensionless.
constexpr double kS1 = 0.0125;  // sec

}  // namespace

std::vector<double> JacksonSpikeGenerator(const std::vector<double> &rates,
                                          double sample_period,
                                          int num_repeats) {
  // Determine the mean of the rate vector. Ignore the negative or zero rates
  // effectively half-rectifying the input signal.
  const uint64_t num_samples = rates.size();
  const double mean_rate = std::accumulate(
      rates.begin(), rates.end(), 0.0, [](double current, double rate) {
        if (rate > 0.0) {
          return current + rate;
        } else {
          return current;
        }
      }) / num_samples;
  const double dt = sample_period;
  const double T = num_samples * dt;  // Total duration of the input (sec).
  uint64_t spike_buffer_size = static_cast<uint64_t>(
      mean_rate * T * num_repeats);
  std::vector<double> spike_times;
  spike_times.reserve(spike_buffer_size);

  // Get a vector of pseudo-random numbers.
  std::random_device rnd_seed;  // Obtain a seed for the random number engine.
  std::mt19937 rnd_gen(rnd_seed());  // Standard Mersenne Twister engine.
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  std::vector<double> rnd_inputs(spike_buffer_size + 1);
  std::generate(rnd_inputs.begin(), rnd_inputs.end(),
                [&dist, &rnd_gen]() { return dist(rnd_gen); });

  // Integer number of discrete time bins within deadtime.
  const uint64_t dead_time_index = std::floor(kDeadTime / dt);

  // Deadtime rounded down to length of an integer number of discrete time bins.
  const double dead_time_rounded = dead_time_index * dt;

  // If y0(t) = kC0 * exp(-t / kS0), then y0(t + dt) = y0(t) * refrac_mult0, the
  // <kC1, kS1> pair is similar.
  const double refrac_mult0 = 1.0 - dt / kS0;
  const double refrac_mult1 = 1.0 - dt / kS1;

  // Calculate effects of a random spike before t=0 on refractoriness and the
  // time-warping sum at t=0.
  uint64_t rnd_buf_index = 0;
  const double end_of_last_dead_time = std::max(
      0.0, std::log(rnd_inputs[rnd_buf_index++]) / rates[0] + kDeadTime);

  // Values of first and second exponenentials in the refractory function.
  double refrac_value0 = kC0 * std::exp(end_of_last_dead_time / kS0);
  double refrac_value1 = kC1 * std::exp(end_of_last_dead_time / kS1);

  // Value of time-warping sum. This is the "integral" of the refractory
  // function normalized by 'dt'.
  double x_sum =
      rates[0] * (-end_of_last_dead_time +
                  kC0 * kS0 * (std::exp(end_of_last_dead_time / kS0) - 1.0) +
                  kC1 * kS1 * (std::exp(end_of_last_dead_time / kS1) - 1.0));

  // Calculate first interspike interval in a homogeneous, unit-rate Poisson
  // process (normalized by 'dt').
  double unit_rate_interval = -std::log(rnd_inputs[rnd_buf_index++]) / dt;

  // NOTE: Both 'unit_rate_interval' and 'x_sum' are divided (or normalized) by
  // 'dt' in order to reduce calculation time. This way we only need to divide
  // by 'dt' once per spike (when calculating 'unit_rate_interval'), instead of
  // multiplying by 'dt' once per time bin (when calculating the new value of
  // 'x_sum').
  double time = dt;
  uint64_t k = 0;
  for (int j = 0; j < num_repeats; ++j) {
    for (; (k < num_samples) && (time < T); ++k, time += dt,
             refrac_value0 *= refrac_mult0,
             refrac_value1 *= refrac_mult1) {
      if (rates[k] <= 0.0) {
        continue;  // Ignore non-positive rates.
      }
      // Update time-warping sum.
      x_sum += (rates[k] * (1.0 - refrac_value0 - refrac_value1));

      // A spike occurs when time-warping sum exceeds interspike "time" in
      // unit-rate process.
      if (x_sum >= unit_rate_interval) {
        spike_times.emplace_back(time);

        // Check if we've reached the capacity for the spike buffer. If so,
        // readjust the buffer and also get more random numbers since these have
        // been used up.
        if (spike_times.size() >= spike_buffer_size) {
          spike_buffer_size += static_cast<uint64_t>(
              std::ceil(mean_rate * ((T - time) + (num_repeats - j - 1) * T)));
          rnd_inputs.resize(spike_buffer_size);
          std::generate(rnd_inputs.begin(), rnd_inputs.end(),
                        [&dist, &rnd_gen]() { return dist(rnd_gen); });
          rnd_buf_index = 0;
        }

        // Next interspike "time" in unit-rate process.
        unit_rate_interval = -std::log(rnd_inputs[rnd_buf_index++]) / dt;
        x_sum = 0.0;

        // Increase index and time to the last time bin in the deadtime, and
        // reset (relative) refractory function.
        k += dead_time_index;
        time += dead_time_rounded;
        refrac_value0 = kC0;
        refrac_value1 = kC1;
      }
    }

    // Reset index and time to begin a new repetion. Don't just set to zero,
    // since new repetition may start within the deadtime following the last
    // spike in the previous repetition.
    time -= T;
    k -= num_samples;
  }

  // Delete spike(s) that occur after the last repetition of the rate function
  // ends.
  while (!spike_times.empty() && (spike_times.back() > T)) {
    spike_times.pop_back();
  }
  std::sort(spike_times.begin(), spike_times.end());
  return spike_times;
}

void JacksonSpikeGenerator(const Eigen::Ref<const Eigen::ArrayXXd> &rates,
                           double sample_period, int num_repeats,
                           Eigen::ArrayXXd *output) {
  const int num_channels = rates.rows();
  const int num_samples = rates.cols();
  *output = Eigen::ArrayXXd::Zero(num_channels, num_samples);
  for (int i = 0; i < num_channels; i++) {
    const auto &channel_vec = rates.row(i);
    const std::vector<double> channel_rates(channel_vec.data(),
                                            channel_vec.data() + num_samples);
    const std::vector<double> &spike_times =
        JacksonSpikeGenerator(channel_rates, sample_period, num_repeats);
    for (const auto spike_time : spike_times) {
      const double spike_sample = spike_time / sample_period;
      (*output)(i, spike_sample) = 1.0;
    }
  }
}

JacksonSpikeGeneratorModel::JacksonSpikeGeneratorModel() :
    AuditoryModel(MODEL_JACKSON_SPIKES, STAGE_AUDITORY_NERVE_SYNAPSE,
                  OUTPUT_SYNAPSE_SPIKE_TIMES, false /* waveform input */) {}

void JacksonSpikeGeneratorModel::Init(const StimulusConfig &stimulus_config) {
  AuditoryModel::Init(stimulus_config);
}

ChannelProperties JacksonSpikeGeneratorModel::GetChannelProperties() const {
  const AuditoryModel *prev_model = previous_model();
  GOOGLE_CHECK(prev_model != nullptr);
  return prev_model->GetChannelProperties();;
}

void JacksonSpikeGeneratorModel::ProcessSegment(
    const Eigen::Ref<const Eigen::ArrayXXd> &input,
    Eigen::ArrayXXd *output) {
  const int sample_rate = stimulus_config().sample_rate();
  GOOGLE_CHECK_LT(0, sample_rate);
  JacksonSpikeGenerator(input,
                        1.0 / sample_rate /* sample period */,
                        1 /* number of repetitions */, output);
}

}  // namespace audition
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
