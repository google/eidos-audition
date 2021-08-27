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

#include "eidos/audition/carfac_model.h"

#include "cpp/agc.h"
#include "cpp/car.h"
#include "cpp/carfac.h"
#include "cpp/ihc.h"
#include "eidos/audition/auditory_model_config.pb.h"
#include "eidos/port/logging.h"

namespace eidos {
namespace audition {
namespace {

// We only support monaural inputs.
constexpr int kNumEars = 1;

// Manufactures CARFAC output configuration given the corresponding
// configuration proto.
CARFACOutput *NewOutput(const CarfacConfig::OutputTypes &output_types) {
  return new CARFACOutput(
      output_types.store_nap(),   // Store Neural Activity Patterns (NAP).
      output_types.store_bm(),    // Store Basilar Membrane (BM) displacements.
      output_types.store_ohc(),   // Store Outer Hair Cells (OHCs).
      output_types.store_agc());  // Store Adaptive Gain Control (AGC).
}

// Default initialization for the Cascade of Asymmetric Resonators (CAR)
// filter parameters. These defaults have been copied from CAR implementation.
void InitDefaults(CarfacConfig::CarParams *car) {
  CARParams def_car;  // Constructed with defaults.
  if (car->velocity_scale() == 0.0) {
    car->set_velocity_scale(def_car.velocity_scale);
  }
  if (car->v_offset() == 0.0) car->set_v_offset(def_car.v_offset);
  if (car->min_zeta() == 0.0) car->set_min_zeta(def_car.min_zeta);
  if (car->max_zeta() == 0.0) car->set_max_zeta(def_car.max_zeta);
  if (car->first_pole_theta() == 0.0) {
    car->set_first_pole_theta(def_car.first_pole_theta);
  }
  if (car->zero_ratio() == 0.0) car->set_zero_ratio(def_car.zero_ratio);
  if (car->high_f_damping_compression() == 0.0) {
    car->set_high_f_damping_compression(def_car.high_f_damping_compression);
  }
  if (car->erb_per_step() == 0.0) car->set_erb_per_step(def_car.erb_per_step);
  if (car->min_pole_hz() == 0.0) car->set_min_pole_hz(def_car.min_pole_hz);
  // The Greenwood map's break frequency in Hertz.
  if (car->erb_break_freq() == 0.0) {
    car->set_erb_break_freq(def_car.erb_break_freq);
  }
  // Glassberg and Moore's high-cf ratio.
  if (car->erb_q() == 0.0) car->set_erb_q(def_car.erb_q);
}

// Default initialization for the IHC filter parameters. These defaults have
// been copied from CAR implementation.
void InitDefaults(CarfacConfig::IhcParams *ihc) {
  IHCParams def_ihc;  // Constructed with default values.
  ihc->set_just_half_wave_rectify(def_ihc.just_half_wave_rectify);  // Force.
  ihc->set_one_capacitor(def_ihc.one_capacitor);  // Force.
  if (ihc->tau_lpf() == 0.0) ihc->set_tau_lpf(def_ihc.tau_lpf);
  if (ihc->tau1_out() == 0.0) ihc->set_tau1_out(def_ihc.tau1_out);
  if (ihc->tau1_in() == 0.0) ihc->set_tau1_in(def_ihc.tau1_in);
  if (ihc->tau2_out() == 0.0) ihc->set_tau2_out(def_ihc.tau2_out);
  if (ihc->tau2_in() == 0.0) ihc->set_tau2_in(def_ihc.tau2_in);
  if (ihc->ac_corner_hz() == 0.0) ihc->set_ac_corner_hz(def_ihc.ac_corner_hz);
}

// Default initialization for the Adaptive Gain Control (AGC) filter parameters.
// These defaults have been copied from CAR implementation.
void InitDefaults(CarfacConfig::AgcParams *agc) {
  AGCParams def_agc;  // Constructed with default values.
  if (agc->num_stages() == 0) agc->set_num_stages(def_agc.num_stages);
  if (agc->agc_stage_gain() == 0.0) {
    agc->set_agc_stage_gain(def_agc.agc_stage_gain);
  }
  if (agc->agc1_scales().empty()) {
    *agc->mutable_agc1_scales() = {
      def_agc.agc1_scales.begin(), def_agc.agc1_scales.end()};
  }
  if (agc->agc2_scales().empty()) {
    *agc->mutable_agc2_scales() = {
      def_agc.agc2_scales.begin(), def_agc.agc2_scales.end()};
  }
  GOOGLE_CHECK_LT(0, agc->agc1_scales().size());
  GOOGLE_CHECK_EQ(agc->agc1_scales().size(), agc->agc2_scales().size());
  if (agc->time_constants().empty()) {
    *agc->mutable_time_constants() = {
      def_agc.time_constants.begin(), def_agc.time_constants.end()};
  }
  GOOGLE_CHECK_EQ(agc->time_constants().size(), agc->agc1_scales().size());
  if (agc->decimation().empty()) {
    *agc->mutable_decimation() = {
      def_agc.decimation.begin(), def_agc.decimation.end()};
  }
  if (agc->agc_mix_coeff() == 0.0) {
    agc->set_agc_mix_coeff(def_agc.agc_mix_coeff);
  }
}

// Default initialization for the configuration.
void InitDefaults(CarfacConfig *config) {
  InitDefaults(config->mutable_car());
  InitDefaults(config->mutable_ihc());
  InitDefaults(config->mutable_agc());

  // If no output signal types were configured, only store neural activity
  // pattern (NAP).
  CarfacConfig::OutputTypes *output_types = config->mutable_output_types();
  if (!output_types->store_nap() && !output_types->store_bm() &&
      !output_types->store_ohc() && !output_types->store_agc()) {
    output_types->set_store_nap(true);
    output_types->set_store_bm(false);
    output_types->set_store_ohc(false);
    output_types->set_store_agc(false);
  }
}

// Initializes CAR parameters from the configuration proto.
void InitFromProto(const CarfacConfig::CarParams &proto, CARParams *car) {
  car->velocity_scale = proto.velocity_scale();
  car->v_offset = proto.v_offset();
  car->min_zeta = proto.min_zeta();
  car->max_zeta = proto.max_zeta();
  car->first_pole_theta = proto.first_pole_theta();
  car->zero_ratio = proto.zero_ratio();
  car->high_f_damping_compression = proto.high_f_damping_compression();
  car->erb_per_step = proto.erb_per_step();
  car->min_pole_hz = proto.min_pole_hz();
  car->erb_break_freq = proto.erb_break_freq();
  car->erb_q = proto.erb_q();
}

// Initializes IHC parameters from the configuration proto.
void InitFromProto(const CarfacConfig::IhcParams &proto, IHCParams *ihc) {
  ihc->just_half_wave_rectify = proto.just_half_wave_rectify();
  ihc->one_capacitor = proto.one_capacitor();
  ihc->tau_lpf = proto.tau_lpf();
  ihc->tau1_out = proto.tau1_out();
  ihc->tau1_in = proto.tau1_in();
  ihc->tau2_out = proto.tau2_out();
  ihc->tau2_in = proto.tau2_in();
  ihc->ac_corner_hz = proto.ac_corner_hz();
}

// Initializes AGC parameters from the configuration proto.
void InitFromProto(const CarfacConfig::AgcParams &proto, AGCParams *agc) {
  agc->num_stages = proto.num_stages();
  agc->agc_stage_gain = proto.agc_stage_gain();
  agc->agc1_scales = {
    proto.agc1_scales().begin(), proto.agc1_scales().end()};
  agc->agc2_scales = {
    proto.agc2_scales().begin(), proto.agc2_scales().end()};
  agc->time_constants = {
    proto.time_constants().begin(), proto.time_constants().end()};
  agc->decimation = {
    proto.decimation().begin(), proto.decimation().end()};
  agc->agc_mix_coeff = proto.agc_mix_coeff();
}

}  // namespace

struct CarfacHandle {
  CarfacHandle(FPType sample_rate, const CARParams &car_params,
               const IHCParams &ihc_params, const AGCParams &agc_params)
      : carfac_(kNumEars, sample_rate, car_params, ihc_params, agc_params) {}
  CarfacHandle() = delete;

  CARFAC carfac_;
};

// Note: The output of the last stage of CARFAC, the neural activity pattern
// (NAP) represents, at least conceptually, the instantaneous firing rates of
// the groups of primary auditory neurons attached to each inner hair cell
// (IHC).
//
// TODO(agutkin): Double-check if NAP can represent firing rates.
CarfacModel::CarfacModel(const CarfacConfig &config) : AuditoryModel(
    MODEL_CARFAC, STAGE_AUDITORY_NERVE_SYNAPSE, OUTPUT_SYNAPSE_FIRING_RATES,
    true /* waveform input */) {
  config_ = config;
  InitDefaults(&config_);
}
CarfacModel::~CarfacModel() = default;

void CarfacModel::Init(const StimulusConfig &stimulus_config) {
  AuditoryModel::Init(stimulus_config);

  // Initialize the values for CARFAC filter parameters.
  CARParams car_params;  // Will use the defaults.
  InitFromProto(config_.car(), &car_params);
  IHCParams ihc_params;  // Will use the defaults.
  InitFromProto(config_.ihc(), &ihc_params);
  AGCParams agc_params;  // Will use the defaults.
  InitFromProto(config_.agc(), &agc_params);

  // Create the filter handle.
  const int sample_rate = stimulus_config.sample_rate();
  GOOGLE_CHECK_LT(0, sample_rate) << "Invalid sampling rate: " << sample_rate;
  impl_.reset(new CarfacHandle(sample_rate, car_params, ihc_params,
                               agc_params));

  // CARFAC configures its own number of filters based on the sampling rate.
  if (stimulus_config.num_channels() > 0) {
    GOOGLE_LOG(WARNING) << "Ignoring override of number of channels";
  }

  // Compute channel properties.
  auto &center_frequencies = channel_properties_.center_frequencies;
  const int num_channels = impl_->carfac_.num_channels();
  const Eigen::ArrayXf &pole_frequencies = impl_->carfac_.pole_frequencies();
  GOOGLE_CHECK_EQ(num_channels, pole_frequencies.size());
  center_frequencies.resize(num_channels);
  Eigen::VectorXd::Map(&center_frequencies[0], num_channels) =
      pole_frequencies.cast<double>();
}

ChannelProperties CarfacModel::GetChannelProperties() const {
  return channel_properties_;
}

void CarfacModel::Reset() {
  GOOGLE_CHECK(impl_ != nullptr);
  impl_->carfac_.Reset();
}

void CarfacModel::ProcessSegment(const Eigen::Ref<const Eigen::ArrayXXd> &input,
                                 Eigen::ArrayXXd *output) {
  GOOGLE_CHECK(impl_ != nullptr);
  std::unique_ptr<CARFACOutput> output_buffer(NewOutput(
      config_.output_types()));
  impl_->carfac_.RunSegment(input.cast<float>(),
                            config_.agc_open_loop(),  // AGC open loop.
                            output_buffer.get());
  const std::vector<ArrayXX> *output_signal = nullptr;
  const CarfacConfig::OutputTypes &output_types = config_.output_types();
  if (output_types.store_nap()) {
    output_signal = &output_buffer->nap();
  } else if (output_types.store_bm()) {
    output_signal = &output_buffer->bm();
  } else if (output_types.store_ohc()) {
    output_signal = &output_buffer->ohc();
  } else {
    output_signal = &output_buffer->agc();
  }
  GOOGLE_CHECK_EQ(1, output_signal->size());  // Single ear.
  *output = (*output_signal)[0].cast<double>();
}

}  // namespace audition
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
