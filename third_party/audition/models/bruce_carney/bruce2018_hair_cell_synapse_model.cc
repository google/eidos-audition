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

#include "third_party/audition/models/bruce_carney/bruce2018_hair_cell_synapse_model.h"

#include <cstdint>
#include <cmath>
#include <vector>

#include "eidos/audition/auditory_model_config.pb.h"
#include "eidos/stubs/logging.h"
#include "src/resample.h"  // See @com_github_resample.
#include "third_party/audition/models/bruce_carney/fractional_gaussian_noise.h"

namespace eidos {
namespace audition {
namespace {

// Sampling frequency used for the synapse model for analysis.
constexpr double kAnalysisSampleRate = 10E3;  // 10 kHz

// Sample period (aka bin width) for the model used in the analysis.
constexpr double kAnalysisSamplePeriod = 1.0 / kAnalysisSampleRate;

// Power-Law function parameters.
constexpr double kAlpha1 = 1.5E-6 * 100E3;
constexpr double kBeta1 = 5E-4;
constexpr double kAlpha2 = 1E-2 * 100E3;
constexpr double kBeta2 = 1E-1;

}  // namespace

// For Power-Law adaptation please see:
//
//   - Drew, P.J. and Abbott, L.F. (2006): "Models and properties of power-law
//     adaptation in neural systems.", Journal of neurophysiology, no. 96,
//     pp. 826-833.
//
// TODO(agutkin): Please note, the implementation below is highly likely to be
// buggy (in terms of container manipulation).
void Bruce2018HairCellSynapse(const Eigen::RowArrayXd &ihc_out,
                              double sample_period,
                              double center_freq, double spont_rate,
                              bool exact_power_law, Eigen::RowArrayXd *rates) {
  const uint32_t num_samples = ihc_out.size();
  const uint32_t delay_point = std::floor(7500.0 / (center_freq / 1E3));
  const uint32_t num_output_samples = std::ceil(
      (num_samples + 2 * delay_point) * sample_period * kAnalysisSampleRate);

  // Generate random numbers. Please note, the mean and variance roughly
  // correspond to the particular implementation of fractional Gaussian noise
  // in the original paper ("ffGn.m").
  double sigma = 0.0;
  if (spont_rate < 0.2) {
    sigma = 1.0;
  } else if (spont_rate < 20.0) {
    sigma = 10.0;
  } else {
    sigma = spont_rate / 2.0;
  }
  const Eigen::VectorXd rand_nums = JacksonFractionalGaussianNoise(
      num_output_samples, 0.9 /* Hurst index */,
      0.0 /* mean */, sigma /* standard deviation */);

  // ================================================================
  // Mapping from the IHC output to the input to Power-Law adapation.
  // ================================================================
  const double cf_slope = std::pow(spont_rate, 0.19) * std::pow(10.0, -0.87);
  const double cf_const = 0.1 * std::pow(std::log10(spont_rate), 2) +
      0.56 * std::log10(spont_rate) - 0.84;
  const double cf_sat = std::pow(10.0, (cf_slope * 8965.5 / 1E3 + cf_const));
  const double cf_factor = std::min(cf_sat, std::pow(
      10.0, cf_slope * center_freq / 1E3 + cf_const)) * 2.0;
  const double mult_factor = std::max(2.95 * std::max(
      1.0, 1.5 - spont_rate / 100), 4.3 - 0.2 * center_freq / 1E3);

  std::vector<double> mapping_out(num_samples);
  for (uint32_t i = 0; i < num_samples; ++i) {
    mapping_out[i] = std::pow(10.0, (0.9 * std::log10(std::fabs(
        ihc_out[i]) * cf_factor)) + mult_factor);
    if (ihc_out[i] < 0.0) mapping_out[i] = -mapping_out[i];
  }

  std::vector<double> power_law_input(num_samples + 3 * delay_point);
  const double spont_rate_3 = 3.0 * spont_rate;
  for (uint32_t i = 0; i < delay_point; ++i) {
    power_law_input[i] = mapping_out[0] + spont_rate_3;
  }
  for (uint32_t i = delay_point; i < num_samples + delay_point; ++i) {
    power_law_input[i] = mapping_out[i - delay_point] + spont_rate_3;
  }
  for (uint32_t i = num_samples + delay_point;
       i < num_samples + 3 * delay_point;
       ++i) {
    power_law_input[i] = power_law_input[i - 1] + spont_rate_3;
  }
  mapping_out.clear();

  // ===============================================================
  // Downsample to low frequency <kAnalysisSampleRate> for analysis.
  // ===============================================================
  const uint32_t downsample_factor = static_cast<uint32_t>(std::ceil(1.0 / (
      sample_period * kAnalysisSampleRate)));
  std::vector<double> resampled_ihc;
  resample(1 /* upsample factor */, downsample_factor, power_law_input,
           resampled_ihc);
  GOOGLE_CHECK_LE(resampled_ihc.size(), power_law_input.size())
      << "Downsampling factor of " << downsample_factor
      << " does not result in smaller output";

  // Pad to <num_output_samples> if necessary.
  if (resampled_ihc.size() < num_output_samples) {
    resampled_ihc.insert(resampled_ihc.end(),
                         power_law_input.begin() + resampled_ihc.size(),
                         power_law_input.begin() + num_output_samples);
  }
  power_law_input.clear();

  // ===============================================================
  // Run Power-Law adaptation.
  // ===============================================================
  uint32_t k = 0;
  double I1 = 0.0, I2 = 0.0;
  std::vector<double> sout1(num_output_samples, 0.0);
  std::vector<double> sout2(num_output_samples, 0.0);
  std::vector<double> n1(num_output_samples, 0.0);
  std::vector<double> n2(num_output_samples, 0.0);
  std::vector<double> n3(num_output_samples, 0.0);
  std::vector<double> m1(num_output_samples, 0.0);
  std::vector<double> m2(num_output_samples, 0.0);
  std::vector<double> m3(num_output_samples, 0.0);
  std::vector<double> m4(num_output_samples, 0.0);
  std::vector<double> m5(num_output_samples, 0.0);
  std::vector<double> resampled_rates(num_output_samples, 0.0);
  for (uint32_t i = 0; i < num_output_samples; ++i, ++k) {
    sout1[k] = std::max(0.0, resampled_ihc[i] + rand_nums[i] - kAlpha1 * I1);
    sout2[k] = std::max(0.0, resampled_ihc[i] - kAlpha2 * I2);
    if (exact_power_law) {
      // Exact implementation of the Power-Law.
      I1 = 0; I2 = 0;
      for (uint32_t j = 0; j <= k; ++j) {
        I1 += (sout1[j] * kAnalysisSamplePeriod /
               ((k - j) * kAnalysisSamplePeriod + kBeta1));
        I2 += (sout2[j] * kAnalysisSamplePeriod /
               ((k - j) * kAnalysisSamplePeriod + kBeta2));
      }
    } else {
      // Approximate implementation of the Power-Law.
      if (k == 0) {
        n1[0] = 1.0E-3 * sout2[0];
        n2[0] = n1[0];
        n3[0] = n2[0];
      } else if (k == 1) {
        n1[1] = 1.992127932802320 * n1[0] + 1.0E-3 * (
            sout2[1] - 0.994466986569624 * sout2[0]);
        n2[1] = 1.999195329360981 * n2[0] + n1[1] - 1.997855276593802 * n1[0];
        n3[1] = -0.798261718183851 * n3[0] + n2[1] + 0.798261718184977 * n2[0];
      } else {
        n1[k] = 1.992127932802320 * n1[k - 1] - 0.992140616993846 * n1[k - 2] +
            1.0E-3 * (sout2[k] - 0.994466986569624 * sout2[k - 1] +
                      0.000000000002347 * sout2[k - 2]);
        n2[k] = 1.999195329360981 * n2[k - 1] - 0.999195402928777 * n2[k - 2] +
            n1[k] - 1.997855276593802 * n1[k - 1] +
            0.997855827934345 * n1[k - 2];
        n3[k] = -0.798261718183851 * n3[k - 1] - 0.199131619873480 * n3[k - 2] +
            n2[k] + 0.798261718184977 * n2[k - 1] +
            0.199131619874064 * n2[k - 2];
      }
      I2 = n3[k];

      if (k == 0) {
        m1[0] = 0.2 * sout1[0];
        m2[0] = m1[0];
        m3[0] = m2[0];
        m4[0] = m3[0];
        m5[0] = m4[0];
      } else if (k == 1) {
        m1[1] = 0.491115852967412 * m1[0] + 0.2 * (
            sout1[1] - 0.173492003319319 * sout1[0]);
        m2[1] = 1.084520302502860 * m2[0] + m1[1] - 0.803462163297112 * m1[0];
        m3[1] = 1.588427084535629 * m3[0] + m2[1] - 1.416084732997016 * m2[0];
        m4[1] = 1.886287488516458 * m4[0] + m3[1] - 1.830362725074550 * m3[0];
        m5[1] = 1.989549282714008 * m5[0] + m4[1] - 1.983165053215032 * m4[0];
      } else {
        m1[k] = 0.491115852967412 * m1[k - 1] - 0.055050209956838 * m1[k - 2] +
            0.2 * (sout1[k] - 0.173492003319319 * sout1[k - 1] +
                   0.000000172983796 * sout1[k - 2]);
        m2[k] = 1.084520302502860 * m2[k - 1] - 0.288760329320566 * m2[k - 2] +
            m1[k] - 0.803462163297112 * m1[k - 1] +
            0.154962026341513 * m1[k - 2];
        m3[k] = 1.588427084535629 * m3[k - 1] - 0.628138993662508 * m3[k - 2] +
            m2[k] - 1.416084732997016 * m2[k - 1] +
            0.496615555008723 * m2[k - 2];
        m4[k] = 1.886287488516458 * m4[k - 1] - 0.888972875389923 * m4[k - 2] +
            m3[k] - 1.830362725074550 * m3[k - 1] +
            0.836399964176882 * m3[k - 2];
        m5[k] = 1.989549282714008 * m5[k - 1] - 0.989558985673023 * m5[k - 2] +
            m4[k] - 1.983165053215032 * m4[k - 1] +
            0.983193027347456 * m4[k - 2];
      }
      I1 = m5[k];
    }
    resampled_rates[k] = sout1[k] + sout2[k];
  }
  GOOGLE_CHECK(k != 0);

  sout1.clear();
  sout2.clear();
  resampled_ihc.clear();

  // ===============================================================
  // Upsample to the original sampling rate.
  // ===============================================================
  std::vector<double> tmp_rates((k - 1) * downsample_factor, 0.0);
  for (uint32_t i = 0; i < k - 1; ++i) {
    const double incr = (resampled_rates[i + 1] - resampled_rates[i]) /
                        downsample_factor;
    for (uint32_t j = 0; j < downsample_factor; ++j) {
      tmp_rates[i * downsample_factor + j] = resampled_rates[i] + j * incr;
    }
  }
  rates->resize(num_samples);
  for (uint32_t i = 0; i < num_samples; ++i) {
    (*rates)[i] = tmp_rates[i + delay_point];
  }
}

Bruce2018HairCellSynapseModel::Bruce2018HairCellSynapseModel() :
    AuditoryModel(MODEL_BRUCE_SYNAPSE_2018, STAGE_AUDITORY_NERVE_SYNAPSE,
                  OUTPUT_SYNAPSE_FIRING_RATES, false /* waveform input */) {}

void Bruce2018HairCellSynapseModel::Init(const StimulusConfig &stimulus_config) {
  AuditoryModel::Init(stimulus_config);
}

ChannelProperties Bruce2018HairCellSynapseModel::GetChannelProperties() const {
  const AuditoryModel *model = previous_model();
  GOOGLE_CHECK(model != nullptr);
  return model->GetChannelProperties();
}

namespace {

// TODO(agutkin): Hard-coded fiber parameters that should be channel-specific.
constexpr double kSpontRate = 60.0;  // Spontaneous firing rate (spikes/sec).

}  // namespace

void Bruce2018HairCellSynapseModel::ProcessSegment(
    const Eigen::Ref<const Eigen::ArrayXXd> &input,
    Eigen::ArrayXXd *output) {
  const int sample_rate = stimulus_config().sample_rate();
  GOOGLE_CHECK_LT(0, sample_rate);

  const ChannelProperties &channel_properties = GetChannelProperties();
  const auto &center_frequencies = channel_properties.center_frequencies;
  const int num_channels = center_frequencies.size();
  GOOGLE_CHECK_LT(0, num_channels);
  *output = Eigen::ArrayXXd::Zero(num_channels, input.cols());
  Eigen::RowArrayXd fiber_output;
  for (int i = 0; i < num_channels; ++i) {
    Bruce2018HairCellSynapse(input.row(i), 1.0 / sample_rate,
                             center_frequencies[i],
                             kSpontRate, false /* approximate Power-Law */,
                             &fiber_output);
    output->row(i).swap(fiber_output.row(0));
  }
}

}  // namespace audition
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
