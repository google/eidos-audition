// Copyright 2007 Ian C. Bruce, Yousof Erfani, Muhammad S. A. Zilany
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

#include "third_party/audition/models/bruce_carney/bruce2018_spike_generator.h"

#include <algorithm>
#include <cstdint>
#include <cmath>
#include <random>

#include "eidos/audition/auditory_model_config.pb.h"
#include "eidos/audition/utils.h"
#include "eidos/stubs/logging.h"

namespace eidos {
namespace audition {
namespace {

constexpr int64_t kZero = static_cast<int64_t>(0);

// -----------------------------------------------
// Synaptic release / spike generation parameters.
// -----------------------------------------------

// Number of synaptic release sites.
constexpr int kNumSites = 4;

// Resting value of the mean redocking time.
constexpr double kTimeRedockingRest = 14.0E-3;

// Size of jump in mean redocking time when a redocking event occurs.
constexpr double kTimeRedockingJump = 0.4E-3;

// Time constant for short-term adaptation (in mean redocking time).
constexpr double kTau = 60.0e-3;

// Minimum spontaneous firing rate (spikes/sec).
constexpr double kMinSpontaneousRate = 1E-4;

// Maximum spontaneous firing rate (spikes/sec).
constexpr double kMaxSpontaneousRate = 180;

// Minimum refractory period (sec).
constexpr double kMinRefractoryPeriod = 0.0;

// Maximum refractory period (sec).
constexpr double kMaxRefractoryPeriod = 20E-3;

// Passes the output of the synapse model through the spike generator.
// The rate probabilities in <rates> correspond to the output of the synapse
// model.
int SpikeGenerator(const std::vector<double> &rates, double sample_period,
                   double t_rd_init, double tabs, double trel,
                   double total_mean_rate, int64_t max_num_spikes,
                   std::vector<double> *spike_times,
                   std::vector<double> *trd_vector) {
  const int64_t num_rates = rates.size();

  // Estimating max number of spikes and events (including before zero
  // elements).
  const double mean_inter_events =
      1.0 / total_mean_rate + t_rd_init / kNumSites;

  // The sufficient array size (more than 99.7% of cases) to register event
  // times after zero.
  const int64_t max_num_events =
      std::ceil(static_cast<int64_t>(
          num_rates * sample_period / mean_inter_events +
          3 * std::sqrt(num_rates * sample_period / mean_inter_events))) +
      kNumSites;

  // Max random array size: <kNumSites> elements for <one_site_redock>
  // initialization, <kNumSites> elements for <pre_release_time_bins>
  // initialization, 1 element for <t_ref> initialization,
  // <max_num_spikes> elements for <t_ref> in the loop,
  // <max_num_events> elements one time for redocking, another time for
  // rate intervals. Also, for before zero elements, average add 2 <kNumSites>
  // events (redock and unit_rate) and add <kNumSites> (max) for <t_ref>s: 3
  // <kNumSites> spikes in total.
  const int64_t rand_buf_len = 2 * kNumSites + 1 + max_num_spikes +
                               2 * max_num_events + max_num_spikes +
                               3 * kNumSites;

  const std::vector<double> rand_nums = GetRandomUniformVector(
      rand_buf_len, 0.0, 1.0);
  int64_t rand_buf_index = 0;

  // Initial < redocking time associated with <kNumSites> release sites.
  std::vector<double> one_site_redock(kNumSites);
  // Initial values associated with nsites release sites.
  std::vector<double> pre_release_time_bins(kNumSites);
  for (int i = 0; i < kNumSites; ++i) {
    one_site_redock[i] = -t_rd_init * std::log(rand_nums[rand_buf_index++]);
    pre_release_time_bins[i] = std::max(
        -static_cast<double>(num_rates),
        std::ceil((kNumSites / std::max(rates[0], 0.1) + t_rd_init) *
                  std::log(rand_nums[rand_buf_index++]) / sample_period));
  }

  // Sort the <kNumSites> initial pre-release times and associate the farthest
  // from zero as the site which has also generated a spike.
  std::sort(pre_release_time_bins.begin(), pre_release_time_bins.end());
  std::vector<double> previous_release_times(kNumSites);
  for (int i = 0; i < kNumSites; ++i) {
    previous_release_times[i] = pre_release_time_bins[i] * sample_period;
  }

  // Position of the first spike (this is often a *negative* value).
  int k_init = static_cast<int>(pre_release_time_bins[0]);
  double current_refractory_period =
      k_init * sample_period;  // Refractory regions.
  // Current refractory time.
  double t_ref = tabs - trel * std::log(rand_nums[rand_buf_index++]);
  int64_t k = k_init;  // Can be negative.

  // Set dynamic mean redocking time to initial mean redocking time.
  double previous_redocking_period = t_rd_init;
  double current_redocking_period = previous_redocking_period;
  // Whether to decay the value of current_redocking_period at the end of the
  // time step.
  bool t_rd_decay = true;
  bool rd_first = false;  // Whether a first redocking event has occurred.

  // Find the spike times for all the num_rates.
  int one_site_redock_rounded, elapsed_time_rounded;
  std::vector<double> elapsed_time(kNumSites, 0.0), x_sum(kNumSites, 0.0);
  std::vector<int> unit_rate_interval(kNumSites, 0);
  std::vector<double> current_release_times(kNumSites, 0.0);
  while (k < num_rates) {
    for (int site_no = 0; site_no < kNumSites; ++site_no) {
      if (k > pre_release_time_bins[site_no]) {
        // Redocking times do not necessarily occur exactly at time step value
        // - calculate the number of integer steps for the elapsed time and
        // redocking time.
        one_site_redock_rounded = static_cast<int>(
            std::floor(one_site_redock[site_no] / sample_period));
        elapsed_time_rounded =
            static_cast<int>(std::floor(elapsed_time[site_no] / sample_period));

        // Check if redocking event has occurred.
        if (one_site_redock_rounded == elapsed_time_rounded) {
          current_redocking_period =
              previous_redocking_period + kTimeRedockingJump;
          previous_redocking_period = current_redocking_period;
          t_rd_decay = false;  // Dont' decay the value of redocking period.
          rd_first = true;     // Redocking has occurred.
        }
        // For each site, the code starts from its associated previous release
        // time.
        elapsed_time[site_no] += sample_period;
      }

      // The elapsed time passes the one time redock (the redocking is
      // finished). In this case the synaptic vesicle starts sensing the input
      // for each site, integration starts after the redockinging is finished
      // for the corresponding site.
      if (elapsed_time[site_no] >= one_site_redock[site_no]) {
        // There are <kNumSites> integrals each vesicle senses 1 / <kNumSites>
        // of the whole rate.
        x_sum[site_no] += (rates[std::max(kZero, k)] / kNumSites);
      }

      // Check if release event happened at <site_no>.
      if ((x_sum[site_no] >= unit_rate_interval[site_no]) &&
          (k >= pre_release_time_bins[site_no])) {
        one_site_redock[site_no] =
            -current_redocking_period * std::log(rand_nums[rand_buf_index++]);
        current_release_times[site_no] =
            previous_release_times[site_no] + elapsed_time[site_no];
        elapsed_time[site_no] = 0.0;

        // Check if a spike occured for the current event - release.
        if (current_release_times[site_no] >= current_refractory_period) {
          // Register only non-negative spike times.
          if (current_release_times[site_no] >= 0) {
            spike_times->push_back(current_release_times[site_no]);
          }
          const double trel_k =
              std::min(trel * 100.0 / rates[std::max(kZero, k)], trel);
          t_ref = tabs - trel_k * std::log(rand_nums[rand_buf_index++]);
          current_refractory_period = current_release_times[site_no] + t_ref;
        }

        previous_release_times[site_no] = current_release_times[site_no];
        x_sum[site_no] = 0.0;
        unit_rate_interval[site_no] = static_cast<int>(
            -std::log(rand_nums[rand_buf_index++]) / sample_period);
      }

      // Error handling in case the random buffer is not large enough.
      if ((rand_buf_index + 1) > rand_buf_len) {
        k = num_rates;
        site_no = kNumSites;
      }
    }

    // Decay the adaptive mean redocking time towards the resting value if no
    // redocking events occurred in this time step.
    if (t_rd_decay && rd_first) {
      current_redocking_period =
          previous_redocking_period -
          (sample_period / kTau) *
              (previous_redocking_period - kTimeRedockingRest);
      previous_redocking_period = current_redocking_period;
    } else {
      t_rd_decay = true;
    }

    // Store the value of the adaptive mean redocking time if it is within the
    // simulation output period.
    if ((k >= 0) && (k < num_rates)) {
      trd_vector->push_back(current_redocking_period);
    }
    k++;
  }
  return spike_times->size();
}

}  // namespace

std::vector<double> Bruce2018SpikeGenerator(const std::vector<double> &rates,
                                            double sample_period,
                                            double spont_rate,
                                            double absolute_refractory_period,
                                            double relative_refractory_period) {
  // Sanity checks.
  CheckWithinRange(spont_rate, kMinSpontaneousRate, kMaxSpontaneousRate,
                   "Spontaneous rate");
  CheckWithinRange(absolute_refractory_period, kMinRefractoryPeriod,
                   kMaxRefractoryPeriod, "Absolute refractory period");
  CheckWithinRange(relative_refractory_period, kMinRefractoryPeriod,
                   kMaxRefractoryPeriod, "Relative refractory period");

  // Initial value of the mean redocking time.
  const double time_redocking_init =
      kTimeRedockingRest + 0.02E-3 * spont_rate - kTimeRedockingJump;

  // Compute the overall mean synaptic rate.
  double total_mean_rate = 0.0;
  const int64_t num_rates = rates.size();
  for (int64_t i = 0; i < num_rates; ++i) {
    total_mean_rate += rates[i] / num_rates;
  }

  // Compute mean interspike interval (ISI) and estimate the number of spikes.
  const double mean_isi =
      1.0 / total_mean_rate + time_redocking_init / kNumSites +
      absolute_refractory_period + relative_refractory_period;
  const double signal_length = num_rates * sample_period;  // In seconds.
  const int64_t max_num_spikes = std::ceil(
      signal_length / mean_isi + 3 * std::sqrt(signal_length / mean_isi));

  // Compute the spikes.
  std::vector<double> spike_times, redocking_times;
  spike_times.reserve(max_num_spikes);
  const int num_spikes = SpikeGenerator(
      rates, sample_period, time_redocking_init, absolute_refractory_period,
      relative_refractory_period, total_mean_rate, max_num_spikes, &spike_times,
      &redocking_times);
  GOOGLE_CHECK_LT(0, num_spikes);
  return spike_times;
}

void Bruce2018SpikeGenerator(const Eigen::Ref<const Eigen::ArrayXXd> &rates,
                             double sample_period, double spont_rate,
                             double absolute_refractory_period,
                             double relative_refractory_period,
                             Eigen::ArrayXXd *output) {
  const int num_channels = rates.rows();
  const int num_samples = rates.cols();
  *output = Eigen::ArrayXXd::Zero(num_channels, num_samples);
  for (int i = 0; i < num_channels; i++) {
    const auto &channel_vec = rates.row(i);
    const std::vector<double> channel_rates(channel_vec.data(),
                                            channel_vec.data() + num_samples);
    const std::vector<double> &spike_times = Bruce2018SpikeGenerator(
        channel_rates, sample_period, spont_rate, absolute_refractory_period,
        relative_refractory_period);
    for (const auto spike_time : spike_times) {
      const double spike_sample = spike_time / sample_period;
      (*output)(i, spike_sample) = 1.0;
    }
  }
}

Bruce2018SpikeGeneratorModel::Bruce2018SpikeGeneratorModel() :
    AuditoryModel(MODEL_BRUCE_SPIKES_2018, STAGE_AUDITORY_NERVE_SYNAPSE,
                  OUTPUT_SYNAPSE_SPIKE_TIMES, false /* waveform input */) {}

void Bruce2018SpikeGeneratorModel::Init(const StimulusConfig &stimulus_config) {
  AuditoryModel::Init(stimulus_config);
}

ChannelProperties Bruce2018SpikeGeneratorModel::GetChannelProperties() const {
  const AuditoryModel *model = previous_model();
  GOOGLE_CHECK(model != nullptr);
  return model->GetChannelProperties();
}

namespace {

// TODO(agutkin): Hard-coded fiber parameters that should be channel-specific.
constexpr double kSpontRate = 60.0;  // Spontaneous firing rate (spikes/sec).
constexpr double kAbsoluteRefractoryPeriod = 0.7E-3;  // sec
constexpr double kRelativeRefractoryPeriod = 0.6E-3;  // sec

}  // namespace

void Bruce2018SpikeGeneratorModel::ProcessSegment(
    const Eigen::Ref<const Eigen::ArrayXXd> &input,
    Eigen::ArrayXXd *output) {
  const int sample_rate = stimulus_config().sample_rate();
  GOOGLE_CHECK_LT(0, sample_rate);
  Bruce2018SpikeGenerator(input,
                          1.0 / sample_rate /* sample period */,
                          kSpontRate, kAbsoluteRefractoryPeriod,
                          kRelativeRefractoryPeriod, output);
}

}  // namespace audition
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
