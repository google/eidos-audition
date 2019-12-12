// Copyright 2000, 2010 Lowel P. O'Mard
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

#include "third_party/audition/models/dsam/sumner2002_hair_cell_synapse_model.h"

#include <float.h>
#include <cmath>
#include <random>
#include <vector>

#include "eidos/stubs/logging.h"

namespace eidos {
namespace audition {
namespace {

struct ChannelVars {
  double act_ca;  // Activation of calcium channels.
  double conc_ca;  // Concentration of calcium at the synapse.
  double reservoir_Q; // Quantity of transmitter in the free pool.
  double cleft_C;  // Quantity of transmitter in the cleft.
  double reprocessed_W;  // Quantity of transm. in the reproc. store.
};

struct UniformRandom {
  UniformRandom() : rnd_gen(rnd_seed()), dist() {}

  std::random_device rnd_seed{};  // Random number seed.
  std::mt19937 rnd_gen;  // Standard Mersenne Twister engine.
  std::uniform_real_distribution<double> dist{};
};

// Following parameters correspond to the HSR population of the fibres
// as published in Table II of the Sumner, et. al. (2002) paper (p. 2183).

// Maximum sample period (sampling interval) the model can hendle (sec).
constexpr double kMaxSamplePeriod = 0.00005;

// Calcium reversal potential, E_Ca (Volts).
constexpr double kCaVRev = 0.066;

// Calcium channel Boltzmann function parameter, beta.
constexpr double kBetaCa = 400.0;

// Calcium channel Boltzmann function parameter, gamma.
constexpr double kGammaCa = 130.0;

// Calcium channel transmitter release exponent (power).
constexpr double kPCa = 3.0;

// Transmitter release permeability, Z (unitless gain).
constexpr double kPerm_z = 2e+32;

// Calcium current time constant (sec).
constexpr double kTauCaChan = 1e-4;

// Calcium ion diffusion (accumulation) time constant (sec).
constexpr double kTauConcCa = 1e-4;

// Replenishment rate (units per second).
constexpr double kReplenishRate_y = 10.0;

// Loss rate (units per second).
constexpr double kLossRate_l = 2580.0;

// Reprocessing rate (units per second).
constexpr double kReprocessRate_x = 66.31;

// Recovery rate (units per second).
constexpr double kRecoveryRate_r = 6580.0;

// HSR population-specific parameters.
namespace hsr {

// Maximum calcium conductance (Siemens).
constexpr double kGCaMax = 8e-09;

// Calcium threshold concentration.
constexpr double kPermCa0 = 4.48e-11;

// Max. no. of transmitter packets in free pool.
constexpr int kMaxFreePool_M = 10;

}  // namespace hsr

// Checks whether the sample period <dt> is good enough for the model.
void CheckSamplePeriod(double dt) {
  GOOGLE_CHECK_LE(dt, kMaxSamplePeriod) << "Maximum sample period exceeded";
  GOOGLE_CHECK_GT(1.0, kReplenishRate_y * dt)
      << "Replenishment rate is too high for the sampling interval";
  GOOGLE_CHECK_GT(1.0, kLossRate_l * dt)
      << "Loss rate is too high for the sampling interval";
  GOOGLE_CHECK_GT(1.0, kReprocessRate_x * dt)
      << "Reprocessing rate is too high for the sampling interval";
  GOOGLE_CHECK_GT(1.0, kRecoveryRate_r * dt)
      << "Recovery rate is too high for the sampling interval";
  GOOGLE_CHECK_GT(1.0, dt / kTauCaChan)
      << "Calcium current time is too high for the sampling interval";
  GOOGLE_CHECK_GT(1.0, dt / kTauConcCa)
      << "Calcium ion diffusion is too high for the sampling interval";
}

// Lifted from DSAM. Note: This is different from std:: implementation which
// counts the number of *failures* before the first success.
//
// This routine calculates the number of times that an event occurs with a set
// probability <prob> and number of trials <num_trials>. It uses the geometric
// distribution method, but the probability must be small.
int GeomDistRandom(double prob, int num_trials, UniformRandom *u_random) {
  if ((num_trials <= 0) || (prob < DBL_EPSILON)) {
    return 0;
  }
  const double log_failure = std::log(1.0 - prob);
  int i;
  double sum;
  for (i = 0, sum = 0.0; ; ++i) {
    const double u_val = u_random->dist(u_random->rnd_gen);
    sum += (std::log(u_val) / log_failure);
    if (sum >= num_trials) break;
  }
  if (i > num_trials) i = num_trials;
  return i;
}

// Checks the configurations and sets the default model parameters if required.
void CheckAndInitDefaults(Sumner2002HairCellSynapseConfig *config) {
  if (config->g_ca_max() <= 0.0) {
    config->set_g_ca_max(hsr::kGCaMax);
  }
  if (config->perm_ca_0() <= 0.0) {
    config->set_perm_ca_0(hsr::kPermCa0);
  }
  if (config->max_free_pool_m() <= 0.0) {
    config->set_max_free_pool_m(hsr::kMaxFreePool_M);
  }
}

}  // namespace

void Sumner2002HairCellSynapse(
    const Sumner2002HairCellSynapseConfig &config,
    const Eigen::Ref<const Eigen::ArrayXXd> &input,
    int sample_rate,
    Eigen::ArrayXXd *output) {
  // Sanity checks.
  GOOGLE_CHECK_LT(0.0, sample_rate);
  const double dt = 1.0 / sample_rate;
  CheckSamplePeriod(dt);

  // Fetch configurable parameters.
  const double ca_cond_mode = config.ca_cond_mode();
  const double cleft_replenish_mode = config.cleft_replenish_mode();
  const bool output_spikes = config.output_spikes();
  const double g_ca_max = config.g_ca_max();
  const double perm_ca0 = config.perm_ca_0();
  const int max_free_pool_m = config.max_free_pool_m();

  // ---------------------------------------------------------------------------
  // Initialization.
  // ---------------------------------------------------------------------------

  // Steady state calcium activation.
  double ssact_ca = 1.0 / (1.0 + std::exp(-kGammaCa * input(0, 0)) / kBetaCa);

  // Calcium conductance (*not* current).
  double i_ca = g_ca_max * std::pow(ssact_ca, 3) * (input(0, 0) - kCaVRev);
  if (ca_cond_mode == SUMNER2002_CA_COND_MODE_REVISED) i_ca *= kTauConcCa;

  double spont_perm_k0 = 0.0;
  if (-i_ca > perm_ca0) {
    spont_perm_k0 = kPerm_z * (std::pow(-i_ca, kPCa) - std::pow(perm_ca0, kPCa));
  }
  double spont_cleft_c0 = kReplenishRate_y / kLossRate_l;
  if (cleft_replenish_mode == SUMNER2002_CLEFT_REPLENISH_MODE_ORIGINAL) {
    spont_cleft_c0 = max_free_pool_m * kReplenishRate_y * spont_perm_k0 /
                     (kReplenishRate_y * (kLossRate_l + kRecoveryRate_r) +
                      spont_perm_k0 * kLossRate_l);
  }
  double spont_free_pool_q0;
  if (spont_cleft_c0 > 0.0) {
    if (!output_spikes) {  // Probabilities.
      spont_free_pool_q0 = spont_cleft_c0 * (kLossRate_l + kRecoveryRate_r) /
                           spont_perm_k0;
    } else {  // Spikes.
      spont_free_pool_q0 =
          std::floor((spont_cleft_c0 * (
              kLossRate_l + kRecoveryRate_r) / spont_perm_k0) + 0.5);
    }
  } else {
    spont_free_pool_q0 = max_free_pool_m;
  }
  const double spont_reprocess_w0 =
      spont_cleft_c0 * kRecoveryRate_r / kReprocessRate_x;

  const int num_channels = input.rows();
  std::vector<ChannelVars> channel_vars(num_channels);
  for (int i = 0; i < num_channels; ++i) {
    channel_vars[i].act_ca = ssact_ca;
    channel_vars[i].conc_ca = -i_ca;
    channel_vars[i].reservoir_Q = spont_free_pool_q0;
    channel_vars[i].cleft_C = spont_cleft_c0;
    channel_vars[i].reprocessed_W = spont_reprocess_w0;
  }

  // ---------------------------------------------------------------------------
  // Produce the estimates given the input.
  // ---------------------------------------------------------------------------
  const int num_samples = input.cols();
  *output = Eigen::ArrayXXd::Zero(num_channels, num_samples);

  const double ydt = kReplenishRate_y * dt;
  const double rdt = kRecoveryRate_r * dt;
  const double xdt = kReprocessRate_x * dt;
  const double zdt = kPerm_z * dt;
  const double k0_pow = std::pow(perm_ca0, kPCa);
  const double l_plus_rdt = (kLossRate_l + kRecoveryRate_r) * dt;
  const double dt_over_tauConcCa = dt / kTauConcCa;
  const double dt_over_tauCaChan = dt / kTauCaChan;
  UniformRandom u_random;   // Uniform random number generator.
  for (int i = 0; i < num_channels; ++i) {
    for (int j = 0; j < num_samples; ++j) {
      // ------------------------------------------------
      // Calcium controlled transmitter release function.
      // ------------------------------------------------
      const double v_in = input(i, j);

      // Calcium current.
      double &act_ca = channel_vars[i].act_ca;
      ssact_ca = 1.0 / (1.0 + std::exp(-kGammaCa * v_in) / kBetaCa);
      act_ca += (ssact_ca - act_ca) * dt_over_tauCaChan;
      i_ca = g_ca_max * std::pow(act_ca, 3) * (v_in - kCaVRev);

      // Calcium ion accumulation and diffusion.
      double &conc_ca = channel_vars[i].conc_ca;
      if (ca_cond_mode == SUMNER2002_CA_COND_MODE_ORIGINAL) {
        conc_ca += (-i_ca - conc_ca) * dt_over_tauConcCa;
      } else {
        // Note: In DSAM, for this particular branch an output-specific, rather
        // than input-specific, sample period is used.
        conc_ca += (-i_ca - conc_ca / kTauConcCa) * dt;
      }

      // Power law release function.
      double kdt = 0.0;
      if (conc_ca > perm_ca0) {
        kdt = zdt * (std::pow(conc_ca, kPCa) - k0_pow);
      }

      // --------
      // Synapse.
      // --------
      if (output_spikes) {   // Spike mode.
        double replenish = 0.0;
        if (cleft_replenish_mode == SUMNER2002_CLEFT_REPLENISH_MODE_UNITY) {
          replenish = GeomDistRandom(ydt, 1, &u_random);
        } else {  // SUMNER2002_CLEFT_REPLENISH_MODE_ORIGINAL
          if (channel_vars[i].reservoir_Q < max_free_pool_m) {
            replenish = GeomDistRandom(
                ydt, static_cast<int>(
                    max_free_pool_m - channel_vars[i].reservoir_Q), &u_random);
          }
        }
        const double reuptake_and_lost = l_plus_rdt * channel_vars[i].cleft_C;
        const double reuptake = rdt * channel_vars[i].cleft_C;

        double reprocessed = 0.0;
        if (channel_vars[i].reprocessed_W >= 1.0) {
          reprocessed = GeomDistRandom(xdt, static_cast<int>(std::floor(
              channel_vars[i].reprocessed_W)), &u_random);
        }
        const double ejected = GeomDistRandom(
            kdt, static_cast<int>(channel_vars[i].reservoir_Q), &u_random);
        channel_vars[i].reservoir_Q += (replenish - ejected + reprocessed);
        channel_vars[i].cleft_C += (ejected - reuptake_and_lost);
        channel_vars[i].reprocessed_W += (reuptake - reprocessed);

        if (ejected > 0.0) {
          (*output)(i, j) = ejected;
        }
      } else {  // Probability mode.
        double replenish = 0.0;
        if (cleft_replenish_mode == SUMNER2002_CLEFT_REPLENISH_MODE_UNITY) {
          replenish = ydt;
        } else {
          if (channel_vars[i].reservoir_Q < max_free_pool_m) {
            replenish = ydt * (max_free_pool_m - channel_vars[i].reservoir_Q);
          }
        }

        const double ejected = kdt * channel_vars[i].reservoir_Q;
        const double reuptake_and_lost = l_plus_rdt * channel_vars[i].cleft_C;
        const double reuptake = rdt * channel_vars[i].cleft_C;
        const double reprocessed = xdt * channel_vars[i].reprocessed_W;

        channel_vars[i].reservoir_Q += (replenish - ejected + reprocessed);
        channel_vars[i].cleft_C += (ejected - reuptake_and_lost);
        channel_vars[i].reprocessed_W += (reuptake - reprocessed);

        (*output)(i, j) = ejected;
      }
    }
  }
}

Sumner2002HairCellSynapseModel::Sumner2002HairCellSynapseModel(
    const Sumner2002HairCellSynapseConfig &config) :
    AuditoryModel(MODEL_SUMNER_SYNAPSE_2002, STAGE_AUDITORY_NERVE_SYNAPSE,
                  OUTPUT_SYNAPSE_FIRING_RATES, false /* waveform input */),
    config_(config) {
  CheckAndInitDefaults(&config_);
}

void Sumner2002HairCellSynapseModel::Init(
    const StimulusConfig &stimulus_config) {
  AuditoryModel::Init(stimulus_config);
}

ChannelProperties Sumner2002HairCellSynapseModel::GetChannelProperties() const {
  const AuditoryModel *prev_model = previous_model();
  GOOGLE_CHECK(prev_model != nullptr);
  return prev_model->GetChannelProperties();;
}

void Sumner2002HairCellSynapseModel::ProcessSegment(
    const Eigen::Ref<const Eigen::ArrayXXd> &input,
    Eigen::ArrayXXd *output) {
  const int sample_rate = stimulus_config().sample_rate();
  GOOGLE_CHECK_LT(0, sample_rate);
  Sumner2002HairCellSynapse(config_, input, sample_rate, output);
}

}  // namespace audition
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
