// Copyright 1999,2000 Frank Baumgarte
// Copyright 2019 Google LLC
//
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file or at
// https://opensource.org/licenses/MIT.

#include "third_party/audition/models/baumgarte/baumgarte_model.h"

#include <cmath>
#include <utility>
#include <vector>

#include "eidos/audition/auditory_model_config.pb.h"
#include "eidos/stubs/integral_types.h"
#include "eidos/stubs/logging.h"

namespace eidos {
namespace audition {
namespace {

// Minimum allowed sample rate in Hz.
constexpr int kMinSampleRate = 100000;

// Maximum allowed number of model sections.
constexpr int kDefaultMaxSections = 252;

// Number of laterally coupled sections in one direction.
constexpr int kMaxLatCoupledSections = 8;

// Highest critical band rate covered by model.
constexpr double kMaxBark = 25.0;

// Model resolution in critical bands.
constexpr double kDeltaZ = 0.1;

// Write every kDownsample value to the output. Should be a positive
// number. By default no downsampling is performed.
constexpr int kDownsample = 1;

// \pi
constexpr double kPi = 3.14159265358979323846;
// \pi * 2
constexpr double kPi2 = 6.283185307179586;

// \pi / 2
constexpr double kHalfPi = 1.5707963267948966;

// Global impedance scaling (for matching to the cochlea input impedance).
constexpr double kRScale = 0.0635;

// Default reference frequency for WDFs [Hz].
constexpr double kDefRefWdfFreq = 1000.0;

// Default reference resistance for WDFs [Ohm].
constexpr double kDefRefWdfRes = 1000.0;

// Amplification used in the linear model.
constexpr double kVV = 100.0;

// Reference voltage: 0.1 mV.
constexpr double kVoltageAmplNorm = 0.0001;

// Audio signal is input at the outer ear or at the oval window.
enum InputSelection { kOvalWindow, kOuterEar };

// Characteristics of non-linearity.
enum NonLinearityCharacteristics {
  kLinear,     // No compression.
  kNonLinear1  // Nonlinear compression of OHC feedback.
};

struct InputEarParameters {
  explicit InputEarParameters(int max_sections) : f_res(max_sections) {}
  InputEarParameters() = delete;

  std::vector<double> f_res;
  int n_sect = -1;
  bool feedback = false;
  NonLinearityCharacteristics nonlin = kNonLinear1;
  bool coupling = false;
  int channel = 0;
  InputSelection inp_sel = kOvalWindow;
  int down_sample = kDownsample;
};

// Structures for hardware model element values.
struct HardwareModelValues {
  explicit HardwareModelValues(int max_sections) :
      c(max_sections), l(max_sections), g(max_sections), q(max_sections),
      ce(max_sections), le(max_sections), ge(max_sections) {}
  HardwareModelValues() = delete;

  double c_q_nu = 0.0, r_q_nu = 0.0;
  double r_hel = 0.0, c_hel = 0.0;
  double r_i = 0.0, z = 0.0;
  // g = 1/r, ce - second amplification stage.
  std::vector<double> c, l, g, q, ce, le, ge;
};

// Structures for WDF resistances of ear model.
struct MiddleEarResistances {
  double gv1, gv2, gv3, gv4, rv5, rv6, rv7, rv8, gv9, gv10, gc2;
};

struct InnerEarResistances {
  explicit InnerEarResistances(int max_sections) :
      rzt(max_sections), g2(max_sections), g4(max_sections), g5(max_sections),
      rL(max_sections), rc(max_sections), r8(max_sections), r9(max_sections),
      r1(max_sections), re2(max_sections), re3(max_sections),
      re1(max_sections) {}
  InnerEarResistances() = delete;

  std::vector<double> rzt, g2, g4, g5, rL, rc, r8, r9, r1;
  double rh1 = 0.0, gh2 = 0.0, gh3 = 0.0;
  std::vector<double> re2, re3, re1;
};

struct WdfResistances {
  explicit WdfResistances(int max_sections) : iear(max_sections) {}
  WdfResistances() = delete;

  MiddleEarResistances mear = {};  // Initialize to defaults.
  InnerEarResistances iear;
};

// Whole hardware model.
struct HardwareModel {
  explicit HardwareModel(int max_sections) : values(max_sections),
                                             wdf_res(max_sections) {}
  HardwareModel() = delete;

  // Can be reset after computing WDF coefficients.
  HardwareModelValues values;
  WdfResistances wdf_res;
};

// Lateral coupling values.
//
// Notes:
//   - "g" is same as in HardwareModelValues.
//   - "rzt" is the same as in InnerEarResistances.
struct CouplingValues {
  explicit CouplingValues(int max_sections) :
      g_l(max_sections), g_k(max_sections), g(max_sections), r(max_sections),
      rzt(max_sections), r_grp(max_sections), g_par(max_sections) {
    for (auto &coupling : g_l) {
      coupling.resize(kMaxLatCoupledSections);
    }
  }
  CouplingValues() = delete;

  std::vector<std::vector<double>> g_l;
  std::vector<double> g_k, g, r, rzt, r_grp, g_par;
};

struct MiddleEarCoefficients {
  double yv2, yv4, yv5, yv6, yv9, yv10;
};

struct InnerEarCoefficients {
  explicit InnerEarCoefficients(int max_sections) :
      y1(max_sections), y4(max_sections), y7(max_sections), y8(max_sections),
      ye1(max_sections), ye2(max_sections) {}
  InnerEarCoefficients() = delete;

  std::vector<double> y1, y4, y7, y8, ye1, ye2;
  double yh1 = 0.0, yh3 = 0.0;
};

struct WdfCoefficients {
  WdfCoefficients(int max_sections) : iear(max_sections) {}
  WdfCoefficients() = delete;

  MiddleEarCoefficients mear = {};  // Initialize to defaults.
  InnerEarCoefficients iear;
};

// Structure for the whole model and analysis containing all coefficients for
// the signal computation.
struct Model {
  explicit Model(int max_sections) : coupl_val(max_sections),
                                     wdf_coef(max_sections) {}
  Model() = delete;

  CouplingValues coupl_val;
  WdfCoefficients wdf_coef;
};

// WDF model signals.
struct Waves {
  explicit Waves(int max_sections) :
      b6(max_sections), a5(max_sections), a7(max_sections), b9(max_sections),
      b2(max_sections), a6(max_sections), b7(max_sections), u(max_sections),
      ae1(max_sections), be3(max_sections), be2(max_sections),
      u_ihc(max_sections), a8(max_sections) {}
  Waves() = delete;

  std::vector<double> b6, a5, a7, b9, b2;
  std::vector<double> a6, b7, u, ae1, be3, be2, u_ihc, a8;
  double bh2 = 0.0, bh1 = 0.0, av1 = 0.0, bv1 = 0.0, bv2 = 0.0, av3 = 0.0,
      bv3 = 0.0, bv5 = 0.0, bv6 = 0.0, av8 = 0.0, bv8 = 0.0, av10 = 0.0,
      bv10 = 0.0;
};

struct VoltagesAndCurrents {
  VoltagesAndCurrents(int max_sections) : u_bm(max_sections),
                                          u_ihc(max_sections) {}
  VoltagesAndCurrents() = delete;

  std::vector<double> u_bm, u_ihc;
  double u_of = 0.0, i_of = 0.0;
  double u_tr = 0.0, i_tr = 0.0;
};

struct InputSignalParameters {
  double sample_freq = 0.0;
};

// Master structure holding all the parameters and state.
struct EarConfigAndState {
  explicit EarConfigAndState(int max_sections) :
      hwm(max_sections), model(max_sections), inp_ear_par(max_sections) {}
  EarConfigAndState() = delete;

  HardwareModel hwm;
  Model model;
  InputSignalParameters inp_sig_par = {};  // Initialize to defaults.
  InputEarParameters inp_ear_par;
};

void InitDefaults(const StimulusConfig &stimulus_config,
                  int native_downsample_step,
                  EarConfigAndState *model) {
  auto *ear_param = &model->inp_ear_par;
  // Signal input at outer ear.
  ear_param->inp_sel = kOuterEar;
  // Number of cochlear model sections.
  // The original implementation sets this to <kMaxBark / kDeltaZ + 1>.
  const int max_recommended_channels = static_cast<int>(kMaxBark / kDeltaZ + 1);
  const int num_channels = stimulus_config.num_channels();
  if (num_channels > max_recommended_channels) {
    GOOGLE_LOG(WARNING) << "Number of requested channels (" << num_channels
                        << ") exceeds number of recommended channels ("
                        << max_recommended_channels << ")";
  }
  ear_param->n_sect = num_channels;
  // Active model with OHC gain.
  ear_param->feedback = true;
  // Nonlinear compression of OHC feedback.
  ear_param->nonlin = kNonLinear1;
  // Coupling of neighbouring sections.
  ear_param->coupling = true;
  // Number of input audio channel (0: 1st channel).
  ear_param->channel = 0;
  // Temporal downsampling factor applied before writing output samples
  // (refers to input sample rate).
  ear_param->down_sample = native_downsample_step;
  // Input sample rate.
  const int sample_rate = stimulus_config.sample_rate();
  if (sample_rate < kMinSampleRate) {
    GOOGLE_LOG(WARNING) << "Sample rate of " << sample_rate
                        << " Hz is too low. Upsample to "
                        << kMinSampleRate << " for sane results";
  }
  model->inp_sig_par.sample_freq = stimulus_config.sample_rate();
}

// Computes values of Zwicker's analog model (inner ear).
void InitHardwareModelValues(InputEarParameters *params,
                             HardwareModelValues *values) {
  const int n_sect = params->n_sect;
  auto &f_res = params->f_res;
  GOOGLE_CHECK_LT(n_sect, f_res.size())
      << "Maximum number of sections exceeded!";

  f_res[1] = 5.0 * kDeltaZ / 0.1;
  f_res[2] = 10.0 * kDeltaZ / 0.1;
  for (int nu = 3; nu <= n_sect; ++nu) {
    f_res[nu] =
        f_res[nu - 1] +
        kDeltaZ * (25.0 + 75.0 * std::pow(1.0 + 1.4E-6 * (f_res[nu - 1] *
                                                          f_res[nu - 1]),
                                          0.69));
  }

  values->c_q_nu = 12.5E-9 * kDeltaZ / 0.1;
  values->r_q_nu = 1.0E6 * kDeltaZ / 0.1;

  auto &c = values->c;
  auto &l = values->l;
  auto &q = values->q;
  auto &g = values->g;
  auto &ce = values->ce;
  auto &le = values->le;
  auto &ge = values->ge;

  for (int nu = 1; nu <= n_sect; ++nu) {
    const double z = nu * kDeltaZ;
    c[nu] = (255.0E-9 + 335.0E-9 * z / 24.0) * kDeltaZ / 0.1;
    l[nu] = 1.0 / (c[nu] * kPi2 * kPi2 * f_res[nu] * f_res[nu]);

    if (z <= 4.0) {
      q[nu] = 1.0 + 0.5 * z;
    }
    if ((z > 4.0) && (z <= 8.0)) {
      q[nu] = 2.0 + 0.25 * z;
    }
    if ((z > 8.0) && (z <= 16.0)) {
      q[nu] = 3.0 + 0.125 * z;
    }
    if (z > 16) {
      q[nu] = 5.0;
    }
    g[nu] = 1.0 / (q[nu] * std::sqrt(l[nu] / c[nu]));

    ce[nu] = c[nu];  // IHC - filter.
    le[nu] = 1.0 / (ce[nu] * kPi2 * kPi2 * 0.9 * 0.9 * f_res[nu] * f_res[nu]);
    ge[nu] = 1.0 / (q[nu] * std::sqrt(le[nu] / ce[nu]));
  }
  values->r_hel = 281.1E3;
  values->c_hel = 56.6E-9;
  values->r_i = 274.0;
}

// Computes resistors for the lateral coupling of sections (inner ear).
void InitResistorsForLateral(HardwareModelValues *values,
                             InputEarParameters *params,
                             CouplingValues *c_val) {
  auto &g_k = c_val->g_k;
  auto &g_l = c_val->g_l;
  auto &r_grp = c_val->r_grp;
  auto &g_par = c_val->g_par;
  auto &rzt = c_val->rzt;
  auto &g = c_val->g;
  auto &r = c_val->r;

  const int n_sect = params->n_sect;
  const double f_coupl = 1.0 / (1.55 * 148.5);
  for (int nu = 1; nu <= n_sect; ++nu) {
    r[nu] = 1.0 / g[nu];
    g_k[nu] = f_coupl * g[nu];
    g_k[nu] *= 0.97;

    const double z = nu * kDeltaZ;
    if (z < 4.5) {
      g_k[nu] *= (1.0 + (4.5 - z) * 0.16);
    }
    if ((z > 10.0) && (z < 17)) {
      g_k[nu] *= (1.0 + (z - 10.0) * 0.023);
    }
    if ((z >= 17.0) && (z < 20)) {
      g_k[nu] *= (1.0 + (z - 10.0) * 0.023 + (z - 17.0) * 0.04);
    }
    if (z >= 20) {
      g_k[nu] *=
          (1.0 + (z - 10.0) * 0.023 + (z - 17.0) * 0.04 + (z - 20.0) * 0.06);
    }
    for (int mu = 1; mu <= kMaxLatCoupledSections; mu++) {
      g_l[nu][mu - 1] =
          g_k[nu] * std::exp(-kHalfPi * mu * mu * kDeltaZ * kDeltaZ * 6.25);
    }
    double tmp = 0.0;
    for (int mu = std::min(kMaxLatCoupledSections, nu - 1); mu > 0; mu--) {
      tmp += g_l[nu][mu - 1];
    }
    for (int mu = std::min(params->n_sect - nu, kMaxLatCoupledSections); mu > 0;
         mu--) {
      tmp += g_l[nu][mu - 1];
    }
    r_grp[nu] = 1.0 / (1.0 + rzt[nu] * (g[nu] + g_k[nu] + tmp));
    g_par[nu] = g_k[nu] + tmp;
  }
}

// Computes WDF resistances for outer/middle ear model.
void InitWdfResistancesMiddleEar(double s_freq, WdfResistances *wdf_resist) {
  const double omega_tmp = kPi * kDefRefWdfFreq / s_freq;
  const double tmp = -kPi2 * kDefRefWdfFreq * std::tan(omega_tmp - kHalfPi);

  auto &mear = wdf_resist->mear;
  mear.gv1 = kDefRefWdfRes / (100000.0 * kRScale);    // Source resistance
  mear.gv2 = kDefRefWdfRes / (tmp * 0.05 * kRScale);  // L1
  mear.gv4 = kDefRefWdfRes / (200.0 * kRScale);       // R1
  mear.gv3 = mear.gv1 + mear.gv2 + mear.gv4;
  mear.rv5 = kRScale / (tmp * 100E-9 * kDefRefWdfRes);   // C2
  mear.rv6 = (tmp * 0.02815 * kRScale) / kDefRefWdfRes;  // L2
  mear.rv7 = (100.0 * kRScale) / kDefRefWdfRes;          // R2
  mear.rv8 = mear.rv5 + mear.rv6 + mear.rv7;
  mear.gv9 = kDefRefWdfRes / (500.0 * kRScale);  // R3
  mear.gv10 = 1.0 / mear.rv8 + mear.gv9;
  mear.gc2 = 1.0 / (1.0 / mear.gv3 + 1.0 / mear.gv10);
}

// Computes WDF resistances for cochlear model (inner ear).
//
// The cot(x) below refers to cotangent of x:
//   cot(x) = 1 / tan(x) = cos(x) / sin(x).
// Since
//   cos(x) = sin(pi/2 - x), sin(x) = cos(pi/2 - x) and tan(-x) = -tan(x),
// the above equation can be rewritten as
//   cot(x) = sin(pi/2 - x) / cos(pi/2 - x) = tan(pi/2 - x) = -tan(x - pi/2).
void InitWdfResistancesInnerEar(double s_freq, HardwareModelValues *values,
                                InputEarParameters *params,
                                WdfResistances *wdf_resist) {
  const auto &f_res = params->f_res;
  const auto &l = values->l;
  const auto &c = values->c;
  const auto &ce = values->ce;
  const auto &le = values->le;
  const auto &ge = values->ge;
  const auto &c_q_nu = values->c_q_nu;
  const auto &r_q_nu = values->r_q_nu;

  const double omega_tmp = kPi / s_freq;
  const int n_sect = params->n_sect;
  auto &i_e = wdf_resist->iear;
  i_e.r1[n_sect] = 1.0 / wdf_resist->mear.gc2;
  for (int nu = n_sect; nu > 0; nu--) {
    const double tmp1 =
        -kPi2 * f_res[nu] *
        std::tan(omega_tmp * f_res[nu] - kHalfPi);  // cot(x)=-tan(x-pi/2)
    i_e.g2[nu] = tmp1 * c_q_nu * kDefRefWdfRes;
    i_e.g4[nu] = kDefRefWdfRes / r_q_nu;
    i_e.g5[nu] = 1.0 / i_e.r1[nu] + i_e.g2[nu] + i_e.g4[nu];
    i_e.rL[nu] = tmp1 * l[nu] / kDefRefWdfRes;
    i_e.rc[nu] = 1.0 / (tmp1 * c[nu] * kDefRefWdfRes);
    i_e.r8[nu] = i_e.rL[nu] * i_e.rc[nu] / (i_e.rL[nu] + i_e.rc[nu]);
    i_e.r9[nu] = i_e.r8[nu] * i_e.rL[nu] / i_e.rc[nu];
    i_e.r1[nu - 1] = 1.0 / i_e.g5[nu] + i_e.r8[nu];

    i_e.rzt[nu] = kDefRefWdfRes * i_e.r8[nu];

    i_e.re2[nu] = 1.0 / (tmp1 * ce[nu] * kDefRefWdfRes);
    i_e.re3[nu] = tmp1 * le[nu] / kDefRefWdfRes;
    i_e.re1[nu] = 1.0 / (ge[nu] * kDefRefWdfRes);
  }

  const double tmp1 =
      -kPi2 * f_res[1] *
      std::tan(omega_tmp * f_res[1] - kHalfPi);  // cot(x)=-tan(x-pi/2)
  wdf_resist->iear.rh1 = i_e.r1[0];
  wdf_resist->iear.gh2 = tmp1 * values->c_hel * kDefRefWdfRes;
  wdf_resist->iear.gh3 = kDefRefWdfRes / values->r_hel;
}

// Computes WDF coefficients.
void InitWdfCoefficients(InputEarParameters *params, WdfResistances *wdf_resist,
                         WdfCoefficients *wdf_coef) {
  const auto &m_e = wdf_resist->mear;
  auto &m_c = wdf_coef->mear;
  m_c.yv2 = m_e.gv2 / m_e.gv3;
  m_c.yv4 = m_e.gv4 / m_e.gv3;
  m_c.yv5 = m_e.rv5 / m_e.rv8;
  m_c.yv6 = m_e.rv6 / m_e.rv8;
  m_c.yv9 = m_e.gv9 / m_e.gv10;
  m_c.yv10 = m_e.gc2 / m_e.gv10;

  auto &i_c = wdf_coef->iear;
  const auto &i_e = wdf_resist->iear;
  const int n_sect = params->n_sect;
  for (int nu = 1; nu <= n_sect; ++nu) {
    i_c.y1[nu] = 1.0 / (i_e.r1[nu] * i_e.g5[nu]);
    i_c.y4[nu] = i_e.g4[nu] / i_e.g5[nu];
    i_c.y7[nu] = i_e.r8[nu] / i_e.r1[nu - 1];
    i_c.y8[nu] = (i_e.r8[nu] - i_e.r9[nu]) / (i_e.r8[nu] + i_e.r9[nu]);

    const double g_tmp =
        1.0 / i_e.re2[nu] + 1.0 / i_e.re3[nu] + 1.0 / i_e.re1[nu];
    i_c.ye1[nu] = 2.0 / (g_tmp * i_e.re1[nu]);
    i_c.ye2[nu] = 2.0 / (g_tmp * i_e.re2[nu]);
  }
  i_c.yh1 = 2.0 / (i_e.rh1 * (1.0 / i_e.rh1 + i_e.gh2 + i_e.gh3));
  i_c.yh3 = 2.0 * i_e.gh3 / (1.0 / i_e.rh1 + i_e.gh2 + i_e.gh3);
}

// Initializes all the model coefficients.
void InitModel(EarConfigAndState *ear_model) {
  auto &inp_sig_par = ear_model->inp_sig_par;
  auto &params = ear_model->inp_ear_par;
  auto &model = ear_model->model;

  const double s_freq = inp_sig_par.sample_freq;
  auto &hwm = ear_model->hwm;
  auto &hwm_val = hwm.values;
  auto &wdf_res = hwm.wdf_res;
  InitHardwareModelValues(&params, &hwm_val);
  InitWdfResistancesMiddleEar(s_freq, &wdf_res);
  InitWdfResistancesInnerEar(s_freq, &hwm_val, &params, &wdf_res);

  auto &coupl_val = model.coupl_val;
  coupl_val.rzt = wdf_res.iear.rzt;
  coupl_val.g = hwm_val.g;

  auto &wdf_coef = model.wdf_coef;
  InitResistorsForLateral(&hwm_val, &params, &coupl_val);
  InitWdfCoefficients(&params, &wdf_res, &wdf_coef);
}

// Implements the WDF algorithm which processes an input audio signal and
// calculates the excitation signals for the inner hair cells along the basilar
// membrane.
void ComputeIHCExcitationWaves(const InputEarParameters &params,
                               const WdfCoefficients &wdf_coef,
                               const CouplingValues &coupl_val, double u_ms,
                               Waves *wn, Waves *wo) {
  // =======================================================
  // Compute nonlinear function (for lateral coupling only):
  // =======================================================
  const int n_sect = params.n_sect;
  std::vector<double> u_nl_vec(n_sect + 1);
  int nu;
  if (params.feedback) {
    switch (params.nonlin) {
      case kLinear: {  // no nonlinearity
        for (nu = 1; nu <= n_sect; ++nu) {
          u_nl_vec[nu] = kVV * wo->u[nu];
        }
        break;
      }
      case kNonLinear1: {  // nonlinearity of middle ear
        if (u_ms < 0) {
          u_ms *= 0.99;
        }
        for (nu = 1; nu <= n_sect; ++nu) {  // nonlinearity of OHC
          const double tmp = wo->u[nu];
          u_nl_vec[nu] =
              100 *
              (tmp / std::sqrt(1.0 + 20.0 * std::fabs(tmp) + 10.0 * tmp * tmp));
        }
        break;
      }
      default: {
        GOOGLE_LOG(FATAL) << "ERROR: unknown nonlinear characteristic";
      }
    }
  }

  // Forward propagating waves.
  const auto &mc = wdf_coef.mear;
  if (params.inp_sel != kOvalWindow) {  // outer ear input
    wn->av1 = u_ms * 600.0;
    wn->av3 = wn->av1 - mc.yv4 * wn->av1 - mc.yv2 * (wn->av1 + wo->bv2);
    wn->av8 = wo->bv6 - wo->bv5;
    wn->av10 = wn->av8 - mc.yv9 * wn->av8;
    wn->b6[n_sect] = -wn->av10 - wn->av3;
  } else {
    wn->b6[n_sect] = u_ms;  // input at oval window
  }
  const auto &ic = wdf_coef.iear;
  for (nu = n_sect; nu > 0; nu--) {
    wn->a5[nu] =
        -(1.0 - ic.y1[nu] - ic.y4[nu]) * wo->b2[nu] - ic.y1[nu] * wn->b6[nu];
    wn->a7[nu] = -wo->b9[nu] - ic.y8[nu] * (wo->a8[nu] + wo->b9[nu]);
    wn->b9[nu] = wn->a7[nu] + wo->a8[nu] + wo->b9[nu];
    wn->b6[nu - 1] = -wn->a5[nu] - wn->a7[nu];
  }
  wn->bh2 = (1.0 - ic.yh3) * wo->bh2 - ic.yh1 * (wo->bh2 - wn->b6[0]);

  // Backward propagating waves.
  wn->a6[0] = wn->bh1 = wn->bh2 + wo->bh2 - wn->b6[0];

  const auto &le = coupl_val;
  const auto &g = coupl_val.g;
  const auto &g_par = coupl_val.g_par;
  const auto &r_grp = coupl_val.r_grp;
  if (!params.feedback) {
    // ==============
    // Passive model:
    // ==============
    for (nu = 1; nu <= n_sect; ++nu) {
      wn->b7[nu] = wn->a7[nu] - ic.y7[nu] * (wn->a6[nu - 1] - wn->b6[nu - 1]);

      wn->u[nu] = wn->b7[nu] / (1.0 + le.rzt[nu] * g[nu]);
      wn->u_ihc[nu] = wn->u[nu];

      wn->a8[nu] = wn->u[nu] + wn->u[nu] - wn->b7[nu];
      wn->b2[nu] = wn->a6[nu - 1] + wn->b7[nu] - wn->a5[nu] - wo->b2[nu];
      wn->a6[nu] = wn->b2[nu] + wo->b2[nu] - wn->b6[nu];
    }
  } else {
    if (!params.coupling) {
      // =========================================
      // Nonlinear model without lateral coupling:
      // =========================================
      for (nu = 1; nu <= n_sect; ++nu) {
        wn->b7[nu] = wn->a7[nu] - ic.y7[nu] * (wn->a6[nu - 1] - wn->b6[nu - 1]);

        wn->u[nu] = (wn->b7[nu] + le.rzt[nu] * le.g_k[nu] * u_nl_vec[nu]) /
                    (1.0 + le.rzt[nu] * (g[nu] + le.g_k[nu]));
        wn->u_ihc[nu] = wn->u[nu];

        wn->a8[nu] = wn->u[nu] + wn->u[nu] - wn->b7[nu];
        wn->b2[nu] = wn->a6[nu - 1] + wn->b7[nu] - wn->a5[nu] - wo->b2[nu];
        wn->a6[nu] = wn->b2[nu] + wo->b2[nu] - wn->b6[nu];
      }
    } else {
      // ======================================
      // Nonlinear model with lateral coupling:
      // ======================================
      for (nu = 1; nu <= n_sect; ++nu) {
        wn->b7[nu] = wn->a7[nu] - ic.y7[nu] * (wn->a6[nu - 1] - wn->b6[nu - 1]);
        double sum_1 = u_nl_vec[nu] * le.g_k[nu];
        for (int mu = std::min(kMaxLatCoupledSections, nu - 1); mu > 0; mu--) {
          sum_1 += u_nl_vec[nu - mu] * le.g_l[nu][mu - 1];
        }
        for (int mu = std::min(params.n_sect - nu, kMaxLatCoupledSections);
             mu > 0; mu--) {
          sum_1 += u_nl_vec[nu + mu] * le.g_l[nu][mu - 1];
        }
        double tmp = sum_1;
        sum_1 *= le.rzt[nu];
        sum_1 += wn->b7[nu];
        wn->u[nu] = sum_1 * r_grp[nu];

        // Second amplification stage.
        wn->ae1[nu] = 20 * (tmp - wn->u[nu] * g_par[nu]) * le.r[nu];

        tmp = wo->be3[nu] + wo->be2[nu];
        wn->be3[nu] = ic.ye1[nu] * (wo->be3[nu] + wn->ae1[nu]) +
                      ic.ye2[nu] * tmp - wo->be3[nu];
        wn->be2[nu] = wn->be3[nu] - tmp;
        wn->u_ihc[nu] = wn->u[nu] + 0.5 * (wn->be2[nu] + wo->be2[nu]);

        wn->a8[nu] = wn->u[nu] + wn->u[nu] - wn->b7[nu];
        wn->b2[nu] = wn->a6[nu - 1] + wn->b7[nu] - wn->a5[nu] - wo->b2[nu];
        wn->a6[nu] = wn->b2[nu] + wo->b2[nu] - wn->b6[nu];
      }
    }
  }

  if (params.inp_sel != kOvalWindow) {
    wn->bv10 = wn->av10 - mc.yv10 * (wn->a6[n_sect] - wn->b6[n_sect]);
    wn->bv3 = -wn->bv10 - wn->a6[n_sect];
    wn->bv8 = wn->bv10 - mc.yv9 * wn->av8;
    wn->bv5 = wo->bv5 - mc.yv5 * (wn->bv8 - wn->av8);
    wn->bv6 = -wo->bv6 - mc.yv6 * (wn->bv8 - wn->av8);
    wn->bv1 = wn->bv3 - mc.yv4 * wn->av1 - mc.yv2 * (wn->av1 + wo->bv2);
    wn->bv2 = wn->bv1 + wn->av1 + wo->bv2;
  }
}

// Transforms wave quantities of the WDF into electrical voltages or currents
// of the electrical prototype filter.
void TransformWavesToVoltagesAndCurrents(const InputEarParameters &params,
                                         const Waves &wn, const Waves &wo,
                                         VoltagesAndCurrents *ui) {
  const int n_sect = params.n_sect;
  for (int nu = 1; nu <= n_sect; ++nu) {
    ui->u_bm[nu] = 0.5 * (wn.a7[nu] + wn.b7[nu]);
    ui->u_ihc[nu] = wn.u_ihc[nu];
  }
}

}  // namespace

// This class hides the gory details involving nested structs that need to be
// parametrized by the number of cochlear sections.
struct BaumgarteHandle {
  BaumgarteHandle(int max_sections) : model(max_sections) {}
  BaumgarteHandle() = delete;

  EarConfigAndState model;
};

BaumgarteModel::BaumgarteModel() :
    AuditoryModel(MODEL_BAUMGARTE, STAGE_HAIR_CELLS,
                  OUTPUT_IHC_TRANSMEMBRANE_POTENTIAL,
                  true /* waveform input */) {}
BaumgarteModel::~BaumgarteModel() = default;

void BaumgarteModel::Init(const StimulusConfig &stimulus_config) {
  AuditoryModel::Init(stimulus_config);

  // Figure out the number of cochlear partitions.
  int max_sections = stimulus_config.num_channels();
  if (max_sections == 0) {  // Use defaults.
    max_sections = kDefaultMaxSections;
  } else {
    max_sections++;
  }
  impl_.reset(new BaumgarteHandle(max_sections));
  mutable_stimulus_config()->set_num_channels(max_sections - 1);

  // Initialize the model.
  InitDefaults(this->stimulus_config(), native_downsample_step_, &impl_->model);
  InitModel(&impl_->model);

  // Compute CF information.
  auto &center_frequencies = channel_properties_.center_frequencies;
  const int num_channels = impl_->model.inp_ear_par.n_sect;
  center_frequencies.reserve(num_channels);
  for (int i = 1; i <= num_channels; ++i) {
    center_frequencies.emplace_back(impl_->model.inp_ear_par.f_res[i]);
  }
}

ChannelProperties BaumgarteModel::GetChannelProperties() const {
  return channel_properties_;
}

void BaumgarteModel::ProcessSegment(
    const Eigen::Ref<const Eigen::ArrayXXd> &input, Eigen::ArrayXXd *output) {
  GOOGLE_CHECK(impl_ != nullptr);

  auto &impl = impl_->model;
  const auto &params = impl.inp_ear_par;
  const auto &model = impl.model;
  const auto &wdf_coef = model.wdf_coef;
  const auto &coupl_val = model.coupl_val;

  const int num_channels = stimulus_config().num_channels();
  const int max_sections = num_channels + 1;
  VoltagesAndCurrents ui(max_sections);
  Waves w_new(max_sections), w_old(max_sections);

  const int down_sample = params.down_sample;
  int down_sample_count = down_sample - 1;
  const int64 num_samples = input.cols();
  const int64 num_output_samples =
      std::ceil(static_cast<double>(num_samples) / down_sample);
  *output = Eigen::ArrayXXd::Zero(num_channels, num_output_samples);
  for (int64 i = 0, n = 0; i < num_samples; ++i) {
    const double u_inp = kVoltageAmplNorm * input(0, i);
    ComputeIHCExcitationWaves(params, wdf_coef, coupl_val, u_inp, &w_new,
                              &w_old);
    down_sample_count++;
    if (down_sample_count == down_sample) {
      TransformWavesToVoltagesAndCurrents(params, w_new, w_old, &ui);
      for (int j = 1; j <= params.n_sect; ++j) {
        (*output)(j - 1, n) = ui.u_ihc[j];
      }
      down_sample_count = 0;
      n++;
    }
    std::swap(w_old, w_new);
  }
}

}  // namespace audition
}  // namespace eidos
