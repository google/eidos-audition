// Copyright 1998 Malcolm Slaney (Interval Research Corporation)
// Copyright 2016 Piotr Majdak
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

#include "third_party/audition/models/amt/meddis1986_hair_cell_synapse_model.h"

#include <cmath>
#include <vector>

#include "eidos/audition/auditory_model_config.pb.h"
#include "eidos/port/logging.h"

namespace eidos {
namespace audition {

// TODO(agutkin): Can be implemented more efficiently using Eigen matrix
// operations.
void Meddis1986HairCellSynapse(const Eigen::Ref<const Eigen::ArrayXXd> &input,
                               int sample_rate, bool subtract_spontaneous_rate,
                               Eigen::ArrayXXd *output) {
  GOOGLE_CHECK_LT(0, sample_rate);

  // Parameters from the original Meddis' papers (see the list references in the
  // header file). Some are rescaled by the sample rate.
  constexpr double A = 5.0;    // Permeability equation.
  constexpr double B = 300.0;  // Permeability equation.
  constexpr double M = 1.0;    // Maximum number of free packets in transmitter
                               // pool.
  const double g = 2000.0;     // Release rate.
  const double y = 5.05;       // Replenishment rate.
  const double l = 2500.0;     // Rate of loss from the cleft.
  const double r = 6580.0;     // Rate of return from the cleft.
  const double x = 66.31;      // Rate of release from reprocessing
                               // to free transmitter.

  // This parameter scales the discharge rate. Adjust as necessary.
  // In combination with the gammatone filterbank (ERBFilterBank),
  // h=50000 will produce a steady-state average discharge
  // probability of about 135 spikes/s within the 1kHz channel,
  // for an input consisting of a 1 kHz sinewave at 60 dB SPL
  // (0 dB SPL corresponds to an RMS level of 1.0 at the
  // input of the gammatone filter).
  constexpr double h = 50000.0;

  // Variable parameters:
  //  - The hair cell contains a quantity of the "free transmitter", q(t), which
  //    leaks through a permeable membrane into a synaptic cleft.
  //  - The permeability, k(t), fluctuates as a function of the instantaneous
  //    amplitude of the acoustic stimulation s(t).
  //  - The synaptic cleft contains a fluctuating amount of transmitter
  //    substance c(t). An amount rc(t)dt is continuously returned to the hair
  //    cell and another amount lc(t)dt is continuously being lost.
  //
  // Compute initial estimates corresponding to spantaneous rate state. The
  // spontaneous firing rate is determined from the equilbrium state of the
  // system with no input: the reservoir contents do not change in time.
  const int num_channels = input.rows();
  const double spont_kt = g * A / (A + B);   // spont. k(t)
  const double spont_cleft = M * y * spont_kt / (l * spont_kt + y * (l + r));
  const double spont_free_pool = spont_cleft * (l + r) / spont_kt;
  const double spont_reprocess = spont_cleft * r / x;
  std::vector<double> ct(num_channels), qt(num_channels), wt(num_channels);
  for (int i = 0; i < num_channels; ++i) {
    ct[i] = spont_cleft;
    qt[i] = spont_free_pool;
    wt[i] = spont_reprocess;
  }

  // Process the input signal.
  const int num_samples = input.cols();
  *output = Eigen::ArrayXXd::Zero(num_channels, num_samples);
  const double dt = 1.0 / sample_rate;
  const double gdt = g * dt;
  const double ydt = y * dt;
  const double ldt = l * dt;
  const double rdt = r * dt;
  const double xdt = x * dt;
  for (int i = 0; i < num_channels; ++i) {
    for (int j = 0; j < num_samples; ++j) {
      const double st = std::max(input(i, j) + A, 0.0);  // s(t)
      const double kdt = gdt * st / (st + B);
      const double replenish = std::max(ydt * (M - qt[i]), 0.0);
      const double eject = kdt * qt[i];
      const double loss = ldt * ct[i];
      const double reuptake = rdt * ct[i];
      const double reprocess = xdt * wt[i];

      qt[i] = std::max(qt[i] + replenish - eject + reprocess, 0.0);
      ct[i] = std::max(ct[i] + eject - loss - reuptake, 0.0);
      wt[i] = std::max(wt[i] + reuptake - reprocess, 0.0);

      // Note: The DSAM implementation produces "probability estimate" by using
      // "h * dt" instead of plain "h" here.
      (*output)(i, j) = ct[i] * h;
      if (subtract_spontaneous_rate) {
        (*output)(i, j) = std::max((*output)(i, j) - spont_cleft, 0.0);
      }
    }
  }
}

Meddis1986HairCellSynapseModel::Meddis1986HairCellSynapseModel() :
    AuditoryModel(MODEL_MEDDIS_SYNAPSE_1986, STAGE_AUDITORY_NERVE_SYNAPSE,
                  OUTPUT_SYNAPSE_FIRING_RATES, false /* waveform input */) {}

void Meddis1986HairCellSynapseModel::Init(
    const StimulusConfig &stimulus_config) {
  AuditoryModel::Init(stimulus_config);
}

ChannelProperties Meddis1986HairCellSynapseModel::GetChannelProperties() const {
  const AuditoryModel *prev_model = previous_model();
  GOOGLE_CHECK(prev_model != nullptr);
  return prev_model->GetChannelProperties();;
}

void Meddis1986HairCellSynapseModel::ProcessSegment(
    const Eigen::Ref<const Eigen::ArrayXXd> &input,
    Eigen::ArrayXXd *output) {
  const int sample_rate = stimulus_config().sample_rate();
  GOOGLE_CHECK_LT(0, sample_rate);
  Meddis1986HairCellSynapse(input, sample_rate,
                            false /* subtract_spontaneous_rate */,
                            output);
}

}  // namespace audition
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
