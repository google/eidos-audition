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

// Various APIs for generating the spike times given the AN synapse rate
// probabilities given Jackson's model.
//
// Based on the original code by Scott Jackson, Copyright (c) 2003.
//
// References:
//    - Jackson B.S., Carney L.H. (2005), "The spontaneous-rate histogram of the
//    auditory nerve can be explained by only two or three spontaneous rates
//    and long-range dependence", J. Assoc. Res. Otolaryngol. 6:148-159.
//    - Zhang X., Heinz M.G., Bruce I.C., Carney L.H. (2001), "A
//    phenomenological model for the responses of auditorynerve fibers: I.
//    Nonlinear tuning with compression and suppression.", J. Acoust. Soc. Am.
//    109:648-670.
//
// The original implementation can be found here:
//   - Documentation:
//   https://www.urmc.rochester.edu/MediaLibraries/URMCMedia/labs/carney-lab/documents/articles/Jackson-SGfast-2003.pdf
//   - Code:
//   https://www.urmc.rochester.edu/MediaLibraries/URMCMedia/labs/carney-lab/codes/SGfast_MEX_files_Updated_08_2011.zip
//
// For an excellent overview see the handouts of Prof. David Heeger:
//   http://www.cns.nyu.edu/~david/handouts/poisson.pdf

#ifndef THIRD_PARTY_AUDITION_MODELS_BRUCE_CARNEY_JACKSON_SPIKE_GENERATOR_H_
#define THIRD_PARTY_AUDITION_MODELS_BRUCE_CARNEY_JACKSON_SPIKE_GENERATOR_H_

#include <vector>

#include "Eigen/Core"
#include "eidos/audition/auditory_model.h"
#include "eidos/audition/stimulus_config.pb.h"

namespace eidos {
namespace audition {

class JacksonSpikeGeneratorModel : public AuditoryModel {
 public:
  JacksonSpikeGeneratorModel();
  ~JacksonSpikeGeneratorModel() override = default;

  // Initializes the model from the model-agnostic configuration.
  void Init(const StimulusConfig &stimulus_config) override;

  // Returns the channel properties for the model.
  ChannelProperties GetChannelProperties() const override;

  // Processes an arbitrary number of samples corresponding to firing rate
  // probabilities, with input of size (number of input channels) by (number of
  // samples).
  //
  // The output is an ArrayXXd of size GetNumberOfFrequencyChannels() by
  // input.cols(). It stores the firing times as fire/no fire boolean mask.
  void ProcessSegment(const Eigen::Ref<const Eigen::ArrayXXd> &input,
                      Eigen::ArrayXXd *output) override;
};

// Generates spike times from an inhomogeneous Poisson process with a fixed
// absolute refractory period and a two-exponential relative refractory period.
//
// This implementation has following important properties:
//   - Uses the time transformation method. This requires one pseudo-random
//     number to be generated per spike.
//   - Uses running approximations to the differential equations of which the
//     exponentials in the relative refractory equation are solutions.
//
// Returns spike times in seconds for an inhomogeneous Poisson process with a
// driving rate specified by the vector <rates> at sampling period of
// <sample_period> (aka <dt>) seconds and with refractoriness. Any negative
// elements in <rates> are treated as zeros. Sampling period <dt> is also used
// as the the width of discrete time bins in the algorithm, and spike times
// will be multiples of <dt>. The input <rates> is repeated <num_repeats>
// contiguous times, but the spike times are given relative to the beginning
// of the repetition (not from the beginning of the first repetition). The
// absolute refractory period (i.e. deadtime) is 0.00075 sec long and the
// relative refractory function is:
//
//   (1 - c0 * exp(-(t - td) / s0) - c1 * exp(-(t - td) / s1),
//
// where 't' is the time since the last spike, 'td' is the length of the
// absolute refractory period, c0 = 0.5, s0 = 0.001 sec, c1 = 0.5, and
// s1 = 0.0125 sec.
//
// The number of repetitions parameter <num_repeats> should match the number
// of repetitions of the stimulus as defined in the hair cell/synaptic models.
std::vector<double> JacksonSpikeGenerator(
    const std::vector<double> &rates,  // Rate probabilities (spikes/sec).
    double sample_period,  // Sampling period of the rate function, in seconds.
    int num_repeats = 1);  // Number of times to repeat the inputs in <rates>.

// Computes and returns per-channel spikes given the per-channel firing rates
// produced by the synapse model. The rows of <rates> correspond to channels
// along the cochlear partition, the columns to the firing rates across time.
// The sampling period of the rate function is given by <sample_period>.
//
// The output has the same dimension as inputs. Each row in the output
// corresponds to the binary indicator of whether a given channel fires or not
// at a given time. The single-channel API above is used to generate this
// information.
//
// The number of repetitions parameter <num_repeats> should match the number
// of repetitions of the stimulus as defined in the hair cell/synaptic models.
//
// TODO(agutkin): Sparse matrix for the outputs is more efficient.
void JacksonSpikeGenerator(const Eigen::Ref<const Eigen::ArrayXXd> &rates,
                           double sample_period,  // sec.
                           int num_repeats, Eigen::ArrayXXd *output);

}  // namespace audition
}  // namespace eidos

#endif  // THIRD_PARTY_AUDITION_MODELS_BRUCE_CARNEY_JACKSON_SPIKE_GENERATOR_H_

// Local Variables:
// mode: c++
// End:
