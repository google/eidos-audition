// Copyright 1998 Malcolm Slaney (Interval Research Corporation)
// Copyright 2006 Ning Ma
// Copyright 2015 João Felipe Santos
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

// Gammatone filterbank that can be used to approximate basilar membrane
// displacements in response to audio stimulus.
//
// This version is based on the original Matlab implementation by Malcolm Slaney
// originally described in: M. Slaney (1993) "An Efficient Implementation of the
// Patterson-Holdsworth Auditory Filter Bank", Apple Tech. Rep. #35. The
// original implementation has since been modified by Malcolm to use a 4th order
// filter. In particular, this implementation is heavily based on the work
// of João Felipe Santos who ported the original code to C++. Please see
// README.txt for more information and licensing.

#ifndef THIRD_PARTY_AUDITION_MODELS_GAMMATONE_GAMMATONE_FILTERBANK_MODEL_H_
#define THIRD_PARTY_AUDITION_MODELS_GAMMATONE_GAMMATONE_FILTERBANK_MODEL_H_

#include <memory>

#include "Eigen/Core"
#include "eidos/audition/auditory_model.h"
#include "eidos/audition/gammatone_filterbank_config.pb.h"
#include "eidos/audition/stimulus_config.pb.h"

namespace eidos {
namespace audition {

struct GammatoneHandle;

class GammatoneFilterbankModel : public AuditoryModel {
 public:
  GammatoneFilterbankModel(const GammatoneFilterbankConfig &config);
  GammatoneFilterbankModel() = delete;
  ~GammatoneFilterbankModel() override;

  // Initializes the model from the model-agnostic configuration.
  // Currently we use the internal defaults.
  void Init(const StimulusConfig &stimulus_config) override;

  // Returns the channel properties for the model.
  ChannelProperties GetChannelProperties() const override;

  // Processes an arbitrary number of monaural audio samples, with input of size
  // (number of input channels or "ears", which is 1) by (number of samples).
  // Note that since the memory layout of an Eigen:ArrayXXd is
  // column-major, this means that all of the values for the first
  // time step are stored before all of the values for the second
  // time step, etc.
  //
  // The output is an ArrayXXd of size GetNumberOfFrequencyChannels() by
  // input.cols(). It stores estimates for basilar membrane displacements.
  void ProcessSegment(const Eigen::Ref<const Eigen::ArrayXXd> &input,
                      Eigen::ArrayXXd *output) override;

  // Default number of processing channels.
  static constexpr int kDefaultNumChannels = 251;

  // Minimal characteristic frequency (in Hz) for the analysis.
  static constexpr double kLowestCFrequency = 50.0;

  // Highest characteristic frequency (in Hz) for the analysis.
  static constexpr double kHighestCFrequency = 10000.0;

 private:
  // Model configuration.
  GammatoneFilterbankConfig config_;

  // Channel properties.
  ChannelProperties channel_properties_;

  // Filterbank implementation.
  std::unique_ptr<GammatoneHandle> impl_;
};

}  // namespace audition
}  // namespace eidos

#endif  // THIRD_PARTY_AUDITION_MODELS_GAMMATONE_GAMMATONE_FILTERBANK_MODEL_H_

// Local Variables:
// mode: c++
// End:
