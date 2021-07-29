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

#include "third_party/audition/models/gammatone/gammatone_filterbank_model.h"

#include <cstdint>
#include <cmath>
#include <complex>
#include <tuple>
#include <vector>

#include "eidos/audition/auditory_model_config.pb.h"
#include "eidos/audition/psychoacoustic_scale_utils.h"
#include "eidos/stubs/logging.h"

namespace eidos {
namespace audition {
namespace {

// ERB parameters. Please see:
//   - B.C.J. Moore and B.R. Glasberg (1983), "Suggested formulae for
//     calculating auditory-filter bandwidths and excitation patterns",
//     J. Acoust. Soc. Am. 74: 750-753.
//   - B.R. Glasberg and B.C.J. Moore (1990), "Derivation of auditory filter
//     shapes from notched-noise data", Hearing Research, 47: 103-138.
constexpr double const kEarQ = 9.26449;
constexpr double const kMinBandwidth = 24.7;
constexpr int const kOrder = 1;

// ERB bandwidth correction for the 4th order filter.
constexpr double kBandwidthCorrection = 1.019;

// Number of filters in the filterbank (per-channel).
const int kNumFilters = 4;

// \pi
constexpr double kPi = 3.14159265358979323846;

// Very small number used for clipping the filter coefficients.
constexpr double kTinyEpsilon = 1E-200;

// Coefficients for a biquadratic digital IIR filter.
//
// This filter implements the transfer function:
//
//   H(z) = (a0 + a1 z^-1 + a2 z^-2) / (1 + b1 z^-1 + b2 z^-2)
//
// See https://en.wikipedia.org/wiki/Digital_biquad_filter.
class BiquadFilter {
 public:
  BiquadFilter() = default;
  ~BiquadFilter() = default;

  // Constructs from matrix representation.
  BiquadFilter(const Eigen::Matrix<double, 5, 1> &coeffs) {
    coeffs_ =
        std::make_tuple(coeffs[0], coeffs[1], coeffs[2], coeffs[3], coeffs[4]);
  }

  // Applies the filter to the input signal.
  Eigen::VectorXd Apply(const Eigen::VectorXd &x) {
    const int num_samples = x.size();
    Eigen::VectorXd y(Eigen::VectorXd::Zero(num_samples));
    double xi, yi, z0 = 0, z1 = 0;
    for (int i = 0; i < num_samples; ++i) {
      xi = x[i];
      yi = std::get<0>(coeffs_) * xi + z0;
      z0 = std::get<1>(coeffs_) * xi - std::get<3>(coeffs_) * yi + z1;
      z1 = std::get<2>(coeffs_) * xi - std::get<4>(coeffs_) * yi;
      y[i] = yi;
    }
    return y;
  }

 private:
  // Coefficients: a0, a1, a2, b1, b2.
  std::tuple<double, double, double, double, double> coeffs_ = {};
};

// Computes an array of N frequencies uniformly spaced between high_freq and
// low_freq on an ERB scale.
//
// See:
//   - Moore, B. C. J., and Glasberg, B. R. (1983). "Suggested formulae for
//     calculating auditory-filter bandwidths and excitation patterns", J.
//     Acoust. Soc. Am. 74, 750-753.
//   - M. Slaney (Apple TR #35), "An Efficient Implementation of the
//     Patterson-Holdsworth Cochlear Filter Bank.", 33-34.
Eigen::VectorXd ERBSpaceSlaney(double low_freq, double high_freq,
                               int num_channels) {
  Eigen::VectorXd cf_array(num_channels);
  constexpr double aux = kEarQ * kMinBandwidth;
  for (int i = 1; i <= num_channels; ++i) {
    cf_array[i - 1] =
        -aux +
        std::exp(i * (-std::log(high_freq + aux) + std::log(low_freq + aux)) /
                 num_channels) * (high_freq + aux);
  }
  return cf_array;
}

// Original comments from Malcom Slaney's Matlab implementation:
//
// This function computes the filter coefficients for a bank of Gammatone
// filters.  These filters were defined by Patterson and Holdworth for
// simulating the cochlea.
//
// The result is returned as an array of filter coefficients.  Each row
// of the filter arrays contains the coefficients for four second order
// filters.  The transfer function for these four filters share the same
// denominator (poles) but have different numerators (zeros). All of these
// coefficients are assembled into one vector that the ERBFilterBank can take
// apart to implement the filter.
//
// The filter bank contains "numChannels" channels that extend from half the
// sampling rate to "low_freq".  Alternatively, if the num_channels input
// argument is a vector, then the values of this vector are taken to be the
// center frequency of each desired filter.  (The low_freq argument is ignored
// in this case.)
std::vector<std::vector<BiquadFilter>> MakeERBFiltersSlaney(
    const StimulusConfig &stimulus_config,
    ChannelProperties *channel_properties) {
  const double sample_rate = stimulus_config.sample_rate();
  const double T = 1.0 / sample_rate;
  const double low_freq = stimulus_config.lowest_cf_hz();
  const double high_freq = stimulus_config.highest_cf_hz();
  const int num_channels = stimulus_config.num_channels();
  const Eigen::VectorXd cf = ERBSpaceSlaney(low_freq, high_freq, num_channels);
  GOOGLE_CHECK_EQ(num_channels, cf.size());
  auto &center_frequencies = channel_properties->center_frequencies;
  center_frequencies.resize(num_channels);
  Eigen::VectorXd::Map(&center_frequencies[0], num_channels) = cf;

  double A0 = T;
  double A2 = 0;
  Eigen::VectorXd B(num_channels);
  Eigen::VectorXd B1(num_channels);
  Eigen::VectorXd B2(num_channels);
  Eigen::VectorXd A11(num_channels);
  Eigen::VectorXd A12(num_channels);
  Eigen::VectorXd A13(num_channels);
  Eigen::VectorXd A14(num_channels);
  Eigen::VectorXd gain(num_channels);
  Eigen::VectorXd erb(num_channels);

  std::vector<std::vector<BiquadFilter>> filterbank(num_channels);
  const std::complex<double> i(0, 1);
  std::complex<double> aux1, aux2, aux3, aux4, aux5, aux6;
  for (int k = 0; k < num_channels; ++k) {
    const double f_cos = std::cos(2.0 * cf[k] * kPi * T);
    const double f_sin = std::sin(2.0 * cf[k] * kPi * T);
    const double sqrt_p = std::sqrt(3.0 + std::pow(2.0, 1.5));
    const double sqrt_n = std::sqrt(3.0 - std::pow(2.0, 1.5));
    const auto f_im1 = std::complex<double>(2) * i * cf[k] * kPi * T;
    const auto f_im2 = std::complex<double>(4) * i * cf[k] * kPi * T;
    const auto f_exp_im2 = std::exp(f_im2);

    erb[k] = std::pow(
        std::pow((cf[k] / kEarQ), kOrder) + std::pow(kMinBandwidth, kOrder),
        1.0 / kOrder);
    B[k] = kBandwidthCorrection * 2 * kPi * erb[k];
    const double exp_p = std::exp(B[k] * T);
    const auto f_exp_im1 = std::exp(-B[k] * T + f_im1);

    B1[k] = -2.0 * f_cos / exp_p;
    B2[k] = std::exp(-2.0 * B[k] * T);
    A11[k] =
        -(-B1[k] * T + 2.0 * sqrt_p * T * f_sin / exp_p) / 2.0;
    A12[k] =
        -(-B1[k] * T - 2.0 * sqrt_p * T * f_sin / exp_p) / 2.0;
    A13[k] =
        -(-B1[k] * T + 2.0 * sqrt_n * T * f_sin / exp_p) / 2.0;
    A14[k] =
        -(-B1[k] * T - 2.0 * sqrt_n * T * f_sin / exp_p) / 2.0;

    aux1 = std::complex<double>(-2) * f_exp_im2 * T;
    aux2 = std::complex<double>(2) * f_exp_im1 * T * (f_cos - sqrt_n * f_sin);
    aux3 = std::complex<double>(2) * f_exp_im1 * T * (f_cos + sqrt_n * f_sin);
    aux4 = std::complex<double>(2) * f_exp_im1 * T * (f_cos - sqrt_p * f_sin);
    aux5 = std::complex<double>(2) * f_exp_im1 * T * (f_cos + sqrt_p * f_sin);
    aux6 = std::pow(-2.0 / std::exp(2 * B[k] * T) - std::complex<double>(2) *
                    f_exp_im2 + std::complex<double>(2) *
                    (std::complex<double>(1) + f_exp_im2) / exp_p, 4);
    gain[k] = std::abs(std::abs((aux1 + aux2) * (aux1 + aux3) * (aux1 + aux4) *
                                (aux1 + aux5) / aux6));
    Eigen::Matrix<double, 5, 1> coeffs1, coeffs2, coeffs3, coeffs4;
    coeffs1 << A0 * (1.0 / gain[k]), A11[k] * (1.0 / gain[k]),
        A2 * (1.0 / gain[k]), B1[k], B2[k];
    coeffs2 << A0, A12[k], A2, B1[k], B2[k];
    coeffs3 << A0, A13[k], A2, B1[k], B2[k];
    coeffs4 << A0, A14[k], A2, B1[k], B2[k];

    std::vector<BiquadFilter> filters(kNumFilters);
    filters[0] = BiquadFilter(coeffs1);
    filters[1] = BiquadFilter(coeffs2);
    filters[2] = BiquadFilter(coeffs3);
    filters[3] = BiquadFilter(coeffs4);
    filterbank[k] = filters;
  }
  return filterbank;
}

// This gammatone filter is based on the implementation by Ning Ma from
// University of Sheffield who, in turn, based his implementation on an
// original algorithm from Martin Cooke's Ph.D thesis (Cooke, 1993) using
// the base-band impulse invariant transformation. This implementation is
// highly efficient in that a mathematical rearrangement is used to
// significantly reduce the cost of computing complex exponentials. For
// more detail on this implementation see
//   http://www.dcs.shef.ac.uk/~ning/resources/gammatone/
//
// Note: Martin Cooke's PhD has been reprinted as M. Cooke (1993): "Modelling
// Auditory Processing and Organisation", Cambridge University Press, Series
// "Distinguished Dissertations in Computer Science", August.
void GammatoneCookeAndMa(const StimulusConfig &stimulus_config,
                         const Eigen::VectorXd &input,  // Input signal.
                         ChannelProperties *channel_properties,
                         Eigen::ArrayXXd *output) {
  const int num_samples = input.size();
  const int num_channels = stimulus_config.num_channels();
  GOOGLE_CHECK_EQ(num_samples, output->cols());  // Samples.
  GOOGLE_CHECK_EQ(num_channels, output->rows());  // Channels.
  auto &center_frequencies = channel_properties->center_frequencies;
  center_frequencies.resize(num_channels);

  // Compute characteristic frequencies equally spaced on ERB scale using
  // the canonical procedure.
  const double low_erb = erb::HzToErbRate(stimulus_config.lowest_cf_hz());
  const double high_erb = erb::HzToErbRate(stimulus_config.highest_cf_hz());
  const double step_erb = (high_erb - low_erb) / (num_channels - 1);
  const int sample_rate = stimulus_config.sample_rate();
  const double tpt = (2.0 * kPi) / sample_rate;
  for (int i = 0; i < num_channels; ++i) {
    // ==================
    // Set-up the filter:
    // ==================
    const double cf = erb::ErbRateToHz(low_erb + i * step_erb);
    center_frequencies[i] = cf;
    const double tptbw = tpt * erb::HzToErb(cf) * kBandwidthCorrection;
    const double a = std::exp(-tptbw);
    // Based on integral of impulse response.
    const double gain = (tptbw * tptbw * tptbw * tptbw) / 3.0;

    // Filter coefficients.
    const double a1 = 4.0 * a;
    const double a2 = -6.0 * a * a;
    const double a3 = 4.0 * a * a * a;
    const double a4 = -a * a * a * a;
    const double a5 = a * a;
    double p1r = 0.0, p2r = 0.0, p3r = 0.0, p4r = 0.0, p1i = 0.0, p2i = 0.0,
           p3i = 0.0, p4i = 0.0;

    // ===========================================================
    // exp(a+i*b) = exp(a)*(cos(b)+i*sin(b))
    // q = exp(-i*tpt*cf*t) = cos(tpt*cf*t) + i*(-sin(tpt*cf*t))
    // qcos = cos(tpt*cf*t)
    // qsin = -sin(tpt*cf*t)
    // ===========================================================
    const double coscf = std::cos(tpt * cf);
    const double sincf = std::sin(tpt * cf );
    double qcos = 1, qsin = 0;   /* t=0 & q = exp(-i*tpt*t*cf)*/
    for (int j = 0; j < num_samples; ++j) {
      // Filter part 1 & shift down to d.c.
      double p0r = qcos * input(j) + a1 * p1r + a2 * p2r + a3 * p3r + a4 * p4r;
      if (p0r < kTinyEpsilon) p0r = 0.0;
      double p0i = qsin * input(j) + a1 * p1i + a2 * p2i + a3 * p3i + a4 * p4i;
      if (p0i < kTinyEpsilon) p0i = 0.0;

      // Filter part 2.
      const double u0r = p0r + a1 * p1r + a5 * p2r;
      const double u0i = p0i + a1 * p1i + a5 * p2i;

      // Update filter results.
      p4r = p3r;
      p3r = p2r;
      p2r = p1r;
      p1r = p0r;
      p4i = p3i;
      p3i = p2i;
      p2i = p1i;
      p1i = p0i;

      // Basilar membrane response
      // 1. Shift up in frequency first:
      //   (u0r+i*u0i) * exp(i*tpt*cf*t) = (u0r+i*u0i) *
      //     (qcos + i*(-qsin))
      // 2. Take the real part only:
      //   bm = real(exp(j*wcf*kT).*u) * gain.
      (*output)(i, j) = (u0r * qcos + u0i * qsin) * gain;

      // The basic idea of saving computational load:
      //   cos(a+b) = cos(a)*cos(b) - sin(a)*sin(b)
      //   sin(a+b) = sin(a)*cos(b) + cos(a)*sin(b)
      //   qcos = cos(tpt*cf*t) = cos(tpt*cf + tpt*cf*(t-1))
      //   qsin = -sin(tpt*cf*t) = -sin(tpt*cf + tpt*cf*(t-1))
      const double old_qcos = qcos;
      qcos = coscf * old_qcos  + sincf * qsin;
      qsin = coscf * qsin - sincf * old_qcos;
    }
  }
}

}  // namespace

// Default number of processing channels.
constexpr int GammatoneFilterbankModel::kDefaultNumChannels;

// Minimal characteristic frequency (in Hz) for the analysis.
constexpr double GammatoneFilterbankModel::kLowestCFrequency;

// Highest characteristic frequency (in Hz) for the analysis.
constexpr double GammatoneFilterbankModel::kHighestCFrequency;

struct GammatoneHandle {
  int num_channels;
  std::vector<std::vector<BiquadFilter>> filterbank;
};

GammatoneFilterbankModel::GammatoneFilterbankModel(
    const GammatoneFilterbankConfig &config) : AuditoryModel(
        MODEL_GAMMATONE_SLANEY, STAGE_BASILAR_MEMBRANE,
        OUTPUT_BASILAR_MEMBRANE_DISPLACEMENT,
        true /* waveform input */), config_(config) {}

GammatoneFilterbankModel::~GammatoneFilterbankModel() = default;

void GammatoneFilterbankModel::Init(const StimulusConfig &stimulus_config) {
  AuditoryModel::Init(stimulus_config);
  impl_.reset(new GammatoneHandle);

  // Check the stimulus configuration and update its parameters if these are not
  // set.
  int num_channels = stimulus_config.num_channels();
  if (num_channels <= 0) {   // Use defaults.
    num_channels = kDefaultNumChannels;
  }
  mutable_stimulus_config()->set_num_channels(num_channels);

  float lowest_cf_hz = stimulus_config.lowest_cf_hz();
  if (lowest_cf_hz <= 0) {  // Choose default.
    lowest_cf_hz = kLowestCFrequency;
  }
  mutable_stimulus_config()->set_lowest_cf_hz(lowest_cf_hz);
  GOOGLE_CHECK_LT(0, lowest_cf_hz);

  const GammatoneFilterbankType filter_type = config_.filter_type();
  float highest_cf_hz = stimulus_config.highest_cf_hz();
  if (highest_cf_hz <= 0) {  // Choose default.
    const double sample_rate = stimulus_config.sample_rate();
    GOOGLE_CHECK_LT(0, sample_rate);
    if (filter_type == GAMMATONE_FILTERBANK_SLANEY) {
      // In Slaney's implementation the filters extend from half the
      // sampling rate to "low_freq".
      highest_cf_hz = sample_rate / 2.0;
    } else {
      highest_cf_hz = kHighestCFrequency;
    }
  }
  mutable_stimulus_config()->set_highest_cf_hz(highest_cf_hz);
  GOOGLE_CHECK_LT(lowest_cf_hz, highest_cf_hz);

  // Filter-specific one-shot initialization.
  if (filter_type == GAMMATONE_FILTERBANK_SLANEY) {
    impl_->filterbank = MakeERBFiltersSlaney(this->stimulus_config(),
                                             &channel_properties_);
    GOOGLE_CHECK_EQ(num_channels, impl_->filterbank.size());
    GOOGLE_CHECK_EQ(kNumFilters, impl_->filterbank[0].size());
  } else if (config_.filter_type() != GAMMATONE_FILTERBANK_COOKE_AND_MA) {
    GOOGLE_LOG(FATAL) << "Only Slaney and Ma's filters are implemented.";
  }
}

ChannelProperties GammatoneFilterbankModel::GetChannelProperties() const {
  return channel_properties_;
}

void GammatoneFilterbankModel::ProcessSegment(
    const Eigen::Ref<const Eigen::ArrayXXd> &input, Eigen::ArrayXXd *output) {
  GOOGLE_CHECK(impl_ != nullptr);

  // Apply filterbank.
  const int64_t num_samples = input.cols();
  const int num_channels = stimulus_config().num_channels();
  *output = Eigen::ArrayXXd::Zero(num_channels, num_samples);
  if (config_.filter_type() == GAMMATONE_FILTERBANK_SLANEY) {
    const Eigen::VectorXd &input_d = input.row(0).cast<double>();
    for (int i = 0; i < num_channels; ++i) {
      const Eigen::VectorXd y1 = impl_->filterbank[i][0].Apply(input_d);
      const Eigen::VectorXd y2 = impl_->filterbank[i][1].Apply(y1);
      const Eigen::VectorXd y3 = impl_->filterbank[i][2].Apply(y2);
      const Eigen::VectorXd y4 = impl_->filterbank[i][3].Apply(y3);
      output->row(i) = y4;
    }
  } else {  // Use Ma's implementation of Cooke's filter.
    GammatoneCookeAndMa(stimulus_config(), input.row(0), &channel_properties_,
                        output);
  }
}

}  // namespace audition
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
