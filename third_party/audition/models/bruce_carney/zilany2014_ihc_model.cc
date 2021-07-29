// Copyright 2013 Muhammad S. A. Zilany, Ian C. Bruce, Rasha A. Ibrahim,
//                Paul C. Nelson and Laurel H. Carney
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

#include "third_party/audition/models/bruce_carney/zilany2014_ihc_model.h"

#include <cstdint>
#include <cmath>
#include <complex>
#include <vector>

#include "eidos/audition/auditory_model_config.pb.h"
#include "eidos/audition/psychoacoustic_scale_utils.h"
#include "eidos/stubs/logging.h"

namespace eidos {
namespace audition {
namespace {

// Type of the species.
enum Species {
  SPECIES_CAT = 0,
  SPECIES_HUMAN_SHERA = 1,    // Shera et al. (PNAS 2002).
  SPECIES_HUMAN_GLASBERG = 2  // Glasberg & Moore (Hear. Res. 1990).
};

// \pi * 2
constexpr double kTwoPi = 6.283185307179586;

// Returns a complex number equal to exp(i * theta). This computes
// (0, theta) -> exp^0 * (sin(theta) + i * cos(theta)).
std::complex<double> CompExp(double theta) {
  return std::exp(std::complex<double>(0, theta));
}

// Gets TauMax, TauMin for the tuning filter. The TauMax is determined by the
// bandwidth / Q10 of the tuning filter at low level. The TauMin is determined
// by the gain change between high and low level.
double GetTauWb(double cf, int species, int order, double *taumax,
                double *taumin) {
  // NOTE(agutkin): The gain should be species-dependent, but it's not.
  double gain, Q10 = 0.0;
  gain = 52.0 / 2.0 * (std::tanh(2.2 * std::log10(cf / 0.6e3) + 0.15) + 1.0);
  if (gain > 60.0) gain = 60.0;
  if (gain < 15.0) gain = 15.0;

  // ratio of TauMin/TauMax according to the gain, order.
  const double ratio = std::pow(10, (-gain / (20.0 * order)));
  if (species == Species::SPECIES_CAT) {
    Q10 = std::pow(10, 0.4708 * std::log10(cf / 1E3) + 0.4664);
  } else if (species == Species::SPECIES_HUMAN_SHERA) {
    Q10 = std::pow((cf / 1000), 0.3) * 12.7 * 0.505 + 0.2085;
  } else if (species == Species::SPECIES_HUMAN_GLASBERG) {
    Q10 = cf / 24.7 / (4.37 * (cf / 1000) + 1) * 0.505 + 0.2085;
  } else {
    GOOGLE_LOG(FATAL) << "Unknown species: " << species;
  }
  const double bw = cf / Q10;
  taumax[0] = 2.0 / (kTwoPi * bw);
  taumin[0] = taumax[0] * ratio;
  return 0;
}

double GetTauBm(double cf, double taumax, double *bmTaumax, double *bmTaumin,
                double *ratio) {
  // NOTE(agutkin): The gain should be species-dependent, but it's not.
  double gain =
      52.0 / 2.0 * (std::tanh(2.2 * std::log10(cf / 0.6E3) + 0.15) + 1.0);
  if (gain > 60.0) gain = 60.0;
  if (gain < 15.0) gain = 15.0;

  constexpr double bwfactor = 0.7;
  constexpr double factor = 2.5;
  ratio[0] = std::pow(10, (-gain / (20.0 * factor)));
  bmTaumax[0] = taumax / bwfactor;
  bmTaumin[0] = bmTaumax[0] * ratio[0];
  return 0.0;
}

// Passes the signal through the signal-path C1 Tenth Order Nonlinear
// Chirp-Gammatone Filter.
//
// NOTE: This function is thread unsafe.
double C1ChirpFilt(double x, double tdres, double cf, int n, double taumax,
                   double rsigma) {
  static double C1gain_norm, C1initphase;
  static double C1input[12][4], C1output[12][4];

  // ======== setup the locations of poles and zeros =======
  const double sigma0 = 1 / taumax;
  const double ipw = 1.01 * cf * kTwoPi - 50;
  const double ipb = 0.2343 * kTwoPi * cf - 1104;
  const double rpa = std::pow(10, std::log10(cf) * 0.9 + 0.55) + 2000;
  const double pzero = std::pow(10, std::log10(cf) * 0.7 + 1.6) + 500;
  double rzero = -pzero;

  // ======================================================

  const int order_of_pole = 10;
  const int half_order_pole = order_of_pole / 2;
  const int order_of_zero = half_order_pole;

  const double fs_bilinear = kTwoPi * cf / std::tan(kTwoPi * cf * tdres / 2);
  const double CF = kTwoPi * cf;

  std::complex<double> p[11];
  double preal, pimg;
  if (n == 0) {
    p[1].real(-sigma0);
    p[1].imag(ipw);
    p[5].real(p[1].real() - rpa);
    p[5].imag(p[1].imag() - ipb);
    p[3].real((p[1].real() + p[5].real()) * 0.5);
    p[3].imag((p[1].imag() + p[5].imag()) * 0.5);
    p[2] = std::conj(p[1]);
    p[4] = std::conj(p[3]);
    p[6] = std::conj(p[5]);
    p[7] = p[1];
    p[8] = p[2];
    p[9] = p[5];
    p[10] = p[6];

    C1initphase = 0.0;
    for (int i = 1; i <= half_order_pole; i++) {
      preal = p[i * 2 - 1].real();
      pimg = p[i * 2 - 1].imag();
      C1initphase = C1initphase + std::atan(CF / (-rzero)) -
                    std::atan((CF - pimg) / (-preal)) -
                    std::atan((CF + pimg) / (-preal));
    }

    // ===================== Initialize C1input & C1output =====================
    for (int i = 1; i <= (half_order_pole + 1); i++) {
      C1input[i][3] = 0;
      C1input[i][2] = 0;
      C1input[i][1] = 0;
      C1output[i][3] = 0;
      C1output[i][2] = 0;
      C1output[i][1] = 0;
    }

    // ===================== normalize the gain =====================

    C1gain_norm = 1.0;
    for (int r = 1; r <= order_of_pole; r++) {
      C1gain_norm = C1gain_norm * (std::pow((CF - p[r].imag()), 2) +
                                   p[r].real() * p[r].real());
    }
  }

  const double norm_gain =
      std::sqrt(C1gain_norm) /
      std::pow(std::sqrt(CF * CF + rzero * rzero), order_of_zero);

  p[1].real(-sigma0 - rsigma);
  if (p[1].real() > 0.0) {
    GOOGLE_LOG(FATAL) << "The system becomes unstable!";
  }

  p[1].imag(ipw);
  p[5].real(p[1].real() - rpa);
  p[5].imag(p[1].imag() - ipb);
  p[3].real((p[1].real() + p[5].real()) * 0.5);
  p[3].imag((p[1].imag() + p[5].imag()) * 0.5);
  p[2] = std::conj(p[1]);
  p[4] = std::conj(p[3]);
  p[6] = std::conj(p[5]);
  p[7] = p[1];
  p[8] = p[2];
  p[9] = p[5];
  p[10] = p[6];

  double phase = 0.0;
  for (int i = 1; i <= half_order_pole; i++) {
    preal = p[i * 2 - 1].real();
    pimg = p[i * 2 - 1].imag();
    phase = phase - std::atan((CF - pimg) / (-preal)) -
            std::atan((CF + pimg) / (-preal));
  }

  rzero = -CF / std::tan((C1initphase - phase) / order_of_zero);
  if (rzero > 0.0) {
    GOOGLE_LOG(FATAL) << "The zeros are in the right-hand plane!";
  }

  // ===================================================================
  // Each iteration of a loop below is for a pair of poles and one zero.
  // ===================================================================

  C1input[1][3] = C1input[1][2];
  C1input[1][2] = C1input[1][1];
  C1input[1][1] = x;

  double dy;
  for (int i = 1; i <= half_order_pole; i++) {
    preal = p[i * 2 - 1].real();
    pimg = p[i * 2 - 1].imag();

    const double temp = std::pow((fs_bilinear - preal), 2) + std::pow(pimg, 2);
    dy = C1input[i][1] * (fs_bilinear - rzero) - 2 * rzero * C1input[i][2] -
         (fs_bilinear + rzero) * C1input[i][3] +
         2 * C1output[i][1] *
             (fs_bilinear * fs_bilinear - preal * preal - pimg * pimg) -
         C1output[i][2] *
             ((fs_bilinear + preal) * (fs_bilinear + preal) + pimg * pimg);
    dy = dy / temp;

    C1input[i + 1][3] = C1output[i][2];
    C1input[i + 1][2] = C1output[i][1];
    C1input[i + 1][1] = dy;
    C1output[i][2] = C1output[i][1];
    C1output[i][1] = dy;
  }
  dy = C1output[half_order_pole][1] * norm_gain;
  // Signal path output is divided by 4 to give correct C1 filter gain.
  const double c1filterout = dy / 4.0;
  return c1filterout;
}

// Parallel path C2 filter: same as the signal-path C1 filter with the OHC
// completely impaired.
//
// NOTE: This function is thread unsafe.
double C2ChirpFilt(double xx, double tdres, double cf, int n, double taumax,
                   double fcohc) {
  static double C2gain_norm, C2initphase;
  static double C2input[12][4], C2output[12][4];
  double preal, pimg;

  // ================ setup the locations of poles and zeros =======

  const double sigma0 = 1 / taumax;
  const double ipw = 1.01 * cf * kTwoPi - 50;
  const double ipb = 0.2343 * kTwoPi * cf - 1104;
  const double rpa = std::pow(10, std::log10(cf) * 0.9 + 0.55) + 2000;
  const double pzero = std::pow(10, std::log10(cf) * 0.7 + 1.6) + 500;
  double rzero = -pzero;

  // ===============================================================

  const int order_of_pole = 10;
  const int half_order_pole = order_of_pole / 2;
  const int order_of_zero = half_order_pole;

  const double fs_bilinear = kTwoPi * cf / std::tan(kTwoPi * cf * tdres / 2);
  const double CF = kTwoPi * cf;

  std::complex<double> p[11];
  if (n == 0) {
    p[1].real(-sigma0);
    p[1].imag(ipw);
    p[5].real(p[1].real() - rpa);
    p[5].imag(p[1].imag() - ipb);
    p[3].real((p[1].real() + p[5].real()) * 0.5);
    p[3].imag((p[1].imag() + p[5].imag()) * 0.5);
    p[2] = std::conj(p[1]);
    p[4] = std::conj(p[3]);
    p[6] = std::conj(p[5]);
    p[7] = p[1];
    p[8] = p[2];
    p[9] = p[5];
    p[10] = p[6];

    C2initphase = 0.0;
    for (int i = 1; i <= half_order_pole; i++) {
      preal = p[i * 2 - 1].real();
      pimg = p[i * 2 - 1].imag();
      C2initphase = C2initphase + std::atan(CF / (-rzero)) -
                    std::atan((CF - pimg) / (-preal)) -
                    std::atan((CF + pimg) / (-preal));
    }

    // ===================== Initialize C2input & C2output =====================
    for (int i = 1; i <= (half_order_pole + 1); i++) {
      C2input[i][3] = 0;
      C2input[i][2] = 0;
      C2input[i][1] = 0;
      C2output[i][3] = 0;
      C2output[i][2] = 0;
      C2output[i][1] = 0;
    }

    // ===================== normalize the gain =====================
    C2gain_norm = 1.0;
    for (int r = 1; r <= order_of_pole; r++) {
      C2gain_norm = C2gain_norm * (std::pow((CF - p[r].imag()), 2) +
                                   p[r].real() * p[r].real());
    }
  }

  const double norm_gain =
      std::sqrt(C2gain_norm) /
      std::pow(std::sqrt(CF * CF + rzero * rzero), order_of_zero);

  p[1].real(-sigma0 * fcohc);
  if (p[1].real() > 0.0) {
    GOOGLE_LOG(FATAL) << "The system becomes unstable!";
  }

  p[1].imag(ipw);
  p[5].real(p[1].real() - rpa);
  p[5].imag(p[1].imag() - ipb);
  p[3].real((p[1].real() + p[5].real()) * 0.5);
  p[3].imag((p[1].imag() + p[5].imag()) * 0.5);
  p[2] = std::conj(p[1]);
  p[4] = std::conj(p[3]);
  p[6] = std::conj(p[5]);
  p[7] = p[1];
  p[8] = p[2];
  p[9] = p[5];
  p[10] = p[6];

  double phase = 0.0;
  for (int i = 1; i <= half_order_pole; i++) {
    preal = p[i * 2 - 1].real();
    pimg = p[i * 2 - 1].imag();
    phase = phase - std::atan((CF - pimg) / (-preal)) -
            std::atan((CF + pimg) / (-preal));
  }

  rzero = -CF / std::tan((C2initphase - phase) / order_of_zero);
  if (rzero > 0.0) {
    GOOGLE_LOG(FATAL) << "The zeros are in the right-hand plane!";
  }

  // ==================================================
  //      time loop begins here
  // ==================================================

  C2input[1][3] = C2input[1][2];
  C2input[1][2] = C2input[1][1];
  C2input[1][1] = xx;

  double dy;
  for (int i = 1; i <= half_order_pole; i++) {
    preal = p[i * 2 - 1].real();
    pimg = p[i * 2 - 1].imag();
    const double temp = std::pow((fs_bilinear - preal), 2) + std::pow(pimg, 2);
    dy = C2input[i][1] * (fs_bilinear - rzero) - 2 * rzero * C2input[i][2] -
         (fs_bilinear + rzero) * C2input[i][3] +
         2 * C2output[i][1] *
             (fs_bilinear * fs_bilinear - preal * preal - pimg * pimg) -
         C2output[i][2] *
             ((fs_bilinear + preal) * (fs_bilinear + preal) + pimg * pimg);
    dy = dy / temp;

    C2input[i + 1][3] = C2output[i][2];
    C2input[i + 1][2] = C2output[i][1];
    C2input[i + 1][1] = dy;

    C2output[i][2] = C2output[i][1];
    C2output[i][1] = dy;
  }

  dy = C2output[half_order_pole][1] * norm_gain;
  const double c2filterout = dy / 4.0;
  return c2filterout;
}

// Passes the signal through the Control path Third Order Nonlinear Gammatone
// Filter.
//
// NOTE: This function is thread unsafe.
double WbGammaTone(double x, double tdres, double centerfreq, int n, double tau,
                   double gain, int order) {
  static double wbphase;
  static std::complex<double> wbgtf[4], wbgtfl[4];

  if (n == 0) {
    wbphase = 0;
    for (int i = 0; i <= order; i++) {
      wbgtfl[i] = CompExp(0) * 0.0;
      wbgtf[i] = CompExp(0) * 0.0;
    }
  }

  const double delta_phase = -kTwoPi * centerfreq * tdres;
  wbphase += delta_phase;

  const double dtmp = tau * 2.0 / tdres;
  const double c1LP = (dtmp - 1) / (dtmp + 1);
  const double c2LP = 1.0 / (dtmp + 1);
  wbgtf[0] = x * CompExp(wbphase);  // frequency shift.

  for (int j = 1; j <= order; j++) {  // IIR Bilinear transformation LPF.
    wbgtf[j] = c2LP * gain * (wbgtf[j - 1] + wbgtfl[j - 1]) + c1LP * wbgtfl[j];
  }
  // frequency shift back up.
  const double out = std::real(CompExp(-wbphase) * wbgtf[order]);
  for (int i = 0; i <= order; i++) {
    wbgtfl[i] = wbgtf[i];
  }
  return out;
}

// Calculates the gain and group delay for the Control path Filter.
double GainGroupDelay(double tdres, double centerfreq, double cf, double tau,
                      int *grdelay) {
  const double tmpcos = std::cos(kTwoPi * (centerfreq - cf) * tdres);
  const double dtmp2 = tau * 2.0 / tdres;
  const double c1LP = (dtmp2 - 1) / (dtmp2 + 1);
  const double c2LP = 1.0 / (dtmp2 + 1);
  const double tmp1 = 1 + c1LP * c1LP - 2 * c1LP * tmpcos;
  const double tmp2 = 2 * c2LP * c2LP * (1 + tmpcos);
  const double wb_gain = std::pow(tmp1 / tmp2, 1.0 / 2.0);
  grdelay[0] = static_cast<int>(
      std::floor((0.5 - (c1LP * c1LP - c1LP * tmpcos) /
                            (1 + c1LP * c1LP - 2 * c1LP * tmpcos))));
  return wb_gain;
}

// Calculates the delay (basilar membrane, synapse, etc. for cat).
double DelayCat(double cf) {
  constexpr double A0 = 3.0;
  constexpr double A1 = 12.5;
  const double x = 11.9 * std::log10(0.80 + cf / 456.0);  // cat mapping.
  const double delay = A0 * std::exp(-x / A1) * 1E-3;
  return delay;
}

// Calculates the delay (basilar membrane, synapse, etc.) for human, based
// on Harte et al. (JASA 2009).
#if 0  // Unused in this version of the model.
double DelayHuman(double cf) {
  constexpr double A = -0.37;
  constexpr double B = 11.09 / 2;
  const double delay = B * std::pow(cf * 1E-3, A) * 1E-3;
  return delay;
}
#endif  // 0

// Gets the output of the OHC Nonlinear Function (Boltzman Function).
//
// The output is normalized to maximum value of 1. The "asym" is the ratio of
// positive Max to negative Max.
double Boltzman(double x, double asym, double s0, double s1, double x1) {
  const double shift = 1.0 / (1.0 + asym);
  const double x0 = s0 * std::log((1.0 / shift - 1) / (1 + std::exp(x1 / s1)));
  const double out1 = 1.0 / (1.0 + std::exp(-(x - x0) / s0) *
                                       (1.0 + std::exp(-(x - x1) / s1))) -
                      shift;
  const double out = out1 / (1 - shift);
  return out;
}

// Gets the output of the OHC Low Pass Filter in the Control path.
//
// NOTE: This function is thread unsafe.
double OhcLowPass(double x, double tdres, double Fc, int n, double gain,
                  int order) {
  static double ohc[4], ohcl[4];
  if (n == 0) {
    for (int i = 0; i < (order + 1); i++) {
      ohc[i] = 0;
      ohcl[i] = 0;
    }
  }

  const double c = 2.0 / tdres;
  const double c1LP = (c - kTwoPi * Fc) / (c + kTwoPi * Fc);
  const double c2LP = kTwoPi * Fc / (kTwoPi * Fc + c);

  ohc[0] = x * gain;
  for (int i = 0; i < order; i++) {
    ohc[i + 1] = c1LP * ohcl[i + 1] + c2LP * (ohc[i] + ohcl[i]);
  }
  for (int j = 0; j <= order; j++) ohcl[j] = ohc[j];
  return ohc[order];
}

// Gets the output of the IHC Low Pass Filter.
//
// NOTE: This function is thread unsafe.
double IhcLowPass(double x, double tdres, double Fc, int n, double gain,
                  int order) {
  static double ihc[8], ihcl[8];
  if (n == 0) {
    for (int i = 0; i < (order + 1); i++) {
      ihc[i] = 0;
      ihcl[i] = 0;
    }
  }

  const double C = 2.0 / tdres;
  const double c1LP = (C - kTwoPi * Fc) / (C + kTwoPi * Fc);
  const double c2LP = kTwoPi * Fc / (kTwoPi * Fc + C);

  ihc[0] = x * gain;
  for (int i = 0; i < order; i++) {
    ihc[i + 1] = c1LP * ihcl[i + 1] + c2LP * (ihc[i] + ihcl[i]);
  }
  for (int j = 0; j <= order; j++) {
    ihcl[j] = ihc[j];
  }
  return ihc[order];
}

// Gets the output of the Control path using Nonlinear Function after OHC.
double NLAfterOhc(double x, double taumin, double taumax, double asym) {
  double minR = 0.05;
  const double R = taumin / taumax;
  if (R < minR) minR = 0.5 * R;

  const double dc = (asym - 1) / (asym + 1.0) / 2.0 - minR;
  const double R1 = R - minR;

  /* This is for new nonlinearity */
  const double s0 = -dc / std::log(R1 / (1 - minR));
  const double x1 = std::fabs(x);
  double out = taumax * (minR + (1.0 - minR) * std::exp(-x1 / s0));
  if (out < taumin) out = taumin;
  if (out > taumax) out = taumax;
  return out;
}

// Gets the output of the IHC Nonlinear Function (Logarithmic Transduction
// Functions).
double NLogarithm(double x, double slope, double asym, double cf) {
  constexpr double corner = 80;
  const double strength = 20.0E6 / std::pow(10, corner / 20);
  double xx = std::log(1.0 + strength * std::fabs(x)) * slope;
  if (x < 0) {
    const double splx = 20 * std::log10(-x / 20E-6);
    const double asym_t = asym - (asym - 1) / (1 + std::exp(splx / 5.0));
    xx = -1 / asym_t * xx;
  }
  return xx;
}

void IhcModel(const Eigen::ArrayXd &px, double cf, int nrep, double tdres,
              int totalstim, double cohc, double cihc, Species species,
              Eigen::RowArrayXd *ihcout) {
  // Allocate dynamic memory for the temporary variables.
  std::vector<double> ihcouttmp(totalstim * nrep, 0.0),
      mey1(totalstim, 0.0),  // middle ear
      mey2(totalstim, 0.0),  // middle ear
      mey3(totalstim, 0.0),  // middle ear
      tmpgain(totalstim, 0.0);

  // Calculate the center frequency for the control-path wideband filter from
  // the location on basilar membrane, based on Greenwood (JASA 1990).
  double bmplace, centerfreq;
  if (species == SPECIES_CAT) {  // Cat frequency shift corresponding to 1.2 mm.
    // Compute location on basilar membrane using CF.
    bmplace = 11.9 * std::log10(0.80 + cf / 456.0);
    // Shift the center frequency.
    centerfreq = 456.0 * (std::pow(10, (bmplace + 1.2) / 11.9) - 0.80);
  } else {  // Human.
    // Human frequency shift corresponding to 1.2 mm.
    // Compute location on basilar membrane using CF.
    bmplace = (35 / 2.1) * std::log10(1.0 + cf / 165.4);
    // Shift the center frequency.
    centerfreq = 165.4 * (std::pow(10, (bmplace + 1.2) / (35 / 2.1)) - 1.0);
  }

  // ====== Parameters for the gain ===========.
  double gain =
      52.0 / 2.0 * (std::tanh(2.2 * std::log10(cf / 0.6e3) + 0.15) + 1.0);
  if (gain > 60.0) gain = 60.0;
  if (gain < 15.0) gain = 15.0;

  // ====== Parameters for the control-path wideband filter =======
  constexpr int bmorder = 3;
  double Taumin[1], Taumax[1];
  GetTauWb(cf, species, bmorder, Taumax, Taumin);

  // ====== Parameters for the signal-path C1 filter ======
  double bmTaumin[1], bmTaumax[1], ratiobm[1];
  GetTauBm(cf, Taumax[0], bmTaumax, bmTaumin, ratiobm);
  const double bmTaubm = cohc * (bmTaumax[0] - bmTaumin[0]) + bmTaumin[0];

  // ====== Parameters for the control-path wideband filter =======
  constexpr int wborder = 3;
  const double TauWBMax = Taumin[0] + 0.2 * (Taumax[0] - Taumin[0]);
  const double TauWBMin = TauWBMax / Taumax[0] * Taumin[0];
  double tauwb = TauWBMax + (bmTaubm - bmTaumax[0]) * (TauWBMax - TauWBMin) /
                                (bmTaumax[0] - bmTaumin[0]);
  int grdelay[1];
  double wbgain = GainGroupDelay(tdres, centerfreq, cf, tauwb, grdelay);
  tmpgain[0] = wbgain;
  double lasttmpgain = wbgain;

  // Nonlinear asymmetry of OHC function and IHC C1 transduction function.
  constexpr double ohcasym = 7.0;
  constexpr double ihcasym = 3.0;

  // Prewarping and related constants for the middle ear.
  constexpr double fp = 1E3;  // Prewarping frequency 1 kHz.
  const double C = kTwoPi * fp / std::tan(kTwoPi / 2 * fp * tdres);
  double m11, m12, m13, m14, m15, m16, m21, m22, m23, m24, m25, m26, m31, m32,
      m33, m34, m35, m36;
  double megainmax;  // max middle-ear gain.
  if (species == SPECIES_CAT) {
    // Cat middle-ear filter - simplified version from Bruce et al. (JASA 2003).
    m11 = C / (C + 693.48);
    m12 = (693.48 - C) / C;
    m13 = 0.0;
    m14 = 1.0;
    m15 = -1.0;
    m16 = 0.0;
    m21 = 1 / (std::pow(C, 2) + 11053 * C + 1.163e8);
    m22 = -2 * std::pow(C, 2) + 2.326e8;
    m23 = std::pow(C, 2) - 11053 * C + 1.163e8;
    m24 = std::pow(C, 2) + 1356.3 * C + 7.4417e8;
    m25 = -2 * std::pow(C, 2) + 14.8834e8;
    m26 = std::pow(C, 2) - 1356.3 * C + 7.4417e8;
    m31 = 1 / (std::pow(C, 2) + 4620 * C + 909059944);
    m32 = -2 * std::pow(C, 2) + 2 * 909059944;
    m33 = std::pow(C, 2) - 4620 * C + 909059944;
    m34 = 5.7585e5 * C + 7.1665e7;
    m35 = 14.333e7;
    m36 = 7.1665e7 - 5.7585e5 * C;
    megainmax = 41.1405;
  } else {  // Human.
    // Human middle-ear filter - based on Pascal et al. (JASA 1998).
    m11 = 1 / (std::pow(C, 2) + 5.9761e3 * C + 2.5255e7);
    m12 = (-2 * std::pow(C, 2) + 2 * 2.5255e7);
    m13 = (std::pow(C, 2) - 5.9761e3 * C + 2.5255e7);
    m14 = (std::pow(C, 2) + 5.6665e3 * C);
    m15 = -2 * std::pow(C, 2);
    m16 = (std::pow(C, 2) - 5.6665e3 * C);
    m21 = 1 / (std::pow(C, 2) + 6.4255e3 * C + 1.3975e8);
    m22 = (-2 * std::pow(C, 2) + 2 * 1.3975e8);
    m23 = (std::pow(C, 2) - 6.4255e3 * C + 1.3975e8);
    m24 = (std::pow(C, 2) + 5.8934e3 * C + 1.7926e8);
    m25 = (-2 * std::pow(C, 2) + 2 * 1.7926e8);
    m26 = (std::pow(C, 2) - 5.8934e3 * C + 1.7926e8);
    m31 = 1 / (std::pow(C, 2) + 2.4891e4 * C + 1.2700e9);
    m32 = (-2 * std::pow(C, 2) + 2 * 1.2700e9);
    m33 = (std::pow(C, 2) - 2.4891e4 * C + 1.2700e9);
    m34 = (3.1137e3 * C + 6.9768e8);
    m35 = 2 * 6.9768e8;
    m36 = (-3.1137e3 * C + 6.9768e8);
    megainmax = 2;
  }

  double meout;
  for (int n = 0; n < totalstim; n++) {
    if (n == 0) {
      // Start of the middle-ear filtering section.
      mey1[0] = m11 * px[0];
      if (species != SPECIES_CAT) mey1[0] = m11 * m14 * px[0];
      mey2[0] = mey1[0] * m24 * m21;
      mey3[0] = mey2[0] * m34 * m31;
      meout = mey3[0] / megainmax;
    } else if (n == 1) {
      mey1[1] = m11 * (-m12 * mey1[0] + px[1] - px[0]);
      if (species != SPECIES_CAT) {
        mey1[1] = m11 * (-m12 * mey1[0] + m14 * px[1] + m15 * px[0]);
      }
      mey2[1] = m21 * (-m22 * mey2[0] + m24 * mey1[1] + m25 * mey1[0]);
      mey3[1] = m31 * (-m32 * mey3[0] + m34 * mey2[1] + m35 * mey2[0]);
      meout = mey3[1] / megainmax;
    } else {
      mey1[n] = m11 * (-m12 * mey1[n - 1] + px[n] - px[n - 1]);
      if (species != SPECIES_CAT) {
        mey1[n] = m11 * (-m12 * mey1[n - 1] - m13 * mey1[n - 2] + m14 * px[n] +
                         m15 * px[n - 1] + m16 * px[n - 2]);
      }
      mey2[n] = m21 * (-m22 * mey2[n - 1] - m23 * mey2[n - 2] + m24 * mey1[n] +
                       m25 * mey1[n - 1] + m26 * mey1[n - 2]);
      mey3[n] = m31 * (-m32 * mey3[n - 1] - m33 * mey3[n - 2] + m34 * mey2[n] +
                       m35 * mey2[n - 1] + m36 * mey2[n - 2]);
      meout = mey3[n] / megainmax;
    }

    // ====================
    // Control-path filter.
    // ====================
    const double wbout1 =
        WbGammaTone(meout, tdres, centerfreq, n, tauwb, wbgain, wborder);
    const double wbout = std::pow((tauwb / TauWBMax), wborder) * wbout1 * 10e3 *
                         std::max(1.0, cf / 5e3);

    // Pass the control signal through OHC Nonlinear Function.
    const double ohcnonlinout = Boltzman(wbout, ohcasym, 12.0, 5.0, 5.0);
    // Lowpass filtering after the OHC nonlinearity.
    const double ohcout = OhcLowPass(ohcnonlinout, tdres, 600, n, 1.0, 2);
    // Nonlinear function after OHC low-pass filter.
    const double tmptauc1 =
        NLAfterOhc(ohcout, bmTaumin[0], bmTaumax[0], ohcasym);
    // Time-constant for the signal-path C1 filter.
    const double tauc1 = cohc * (tmptauc1 - bmTaumin[0]) + bmTaumin[0];
    // Shift of the location of poles of the C1 filter from the initial
    // positions.
    const double rsigma = 1.0 / tauc1 - 1.0 / bmTaumax[0];
    if (1.0 / tauc1 < 0.0) {
      GOOGLE_LOG(FATAL) << "The poles are in the right-half plane; "
                        << "system is unstable.";
    }
    tauwb = TauWBMax + (tauc1 - bmTaumax[0]) * (TauWBMax - TauWBMin) /
                           (bmTaumax[0] - bmTaumin[0]);
    const double wb_gain =
        GainGroupDelay(tdres, centerfreq, cf, tauwb, grdelay);
    const int grd = grdelay[0];
    if ((grd + n) < totalstim) tmpgain[grd + n] = wb_gain;
    if (tmpgain[n] == 0) tmpgain[n] = lasttmpgain;
    wbgain = tmpgain[n];
    lasttmpgain = wbgain;

    // ====== Signal-path C1 filter ======.
    const double c1filterouttmp =
        C1ChirpFilt(meout, tdres, cf, n, bmTaumax[0], rsigma);

    // ====== Parallel-path C2 filter ======.
    const double c2filterouttmp =
        C2ChirpFilt(meout, tdres, cf, n, bmTaumax[0], 1.0 / ratiobm[0]);

    // === Run the inner hair cell (IHC) section: NL function and then
    // perform lowpass filtering ===
    const double c1vihctmp =
        NLogarithm(cihc * c1filterouttmp, 0.1, ihcasym, cf);
    const double c2vihctmp = -NLogarithm(
        c2filterouttmp * std::fabs(c2filterouttmp) * cf / 10.0 * cf / 2e3, 0.2,
        1.0, cf);  // C2 transduction output.
    ihcouttmp[n] = IhcLowPass(c1vihctmp + c2vihctmp, tdres, 3000, n, 1.0, 7);
  }

  // Stretch out the IHC output according to nrep (number of repetitions).
  for (int i = 0; i < totalstim * nrep; i++) {
    ihcouttmp[i] = ihcouttmp[static_cast<int>(std::fmod(i, totalstim))];
  }

  // Adjust total path delay to IHC output signal.
  const double delay = DelayCat(cf);  // Human parameters not used in v5.2.
  const int delaypoint =
      std::max(0, static_cast<int>(std::ceil(delay / tdres)));
  for (int i = delaypoint; i < totalstim * nrep; i++) {
    (*ihcout)(i) = ihcouttmp[i - delaypoint];
  }
}

}  // namespace

Zilany2014IhcModel::Zilany2014IhcModel() :
    AuditoryModel(MODEL_ZILANY_IHC_2014, STAGE_HAIR_CELLS,
                  OUTPUT_IHC_TRANSMEMBRANE_POTENTIAL,
                  true /* waveform input */) {}

namespace {

// Default number of cochlear channels.
constexpr int kDefaultNumChannels = 251;

// Outer Hair Cell (OHC) impairment coefficient. Should be between 0
// (degenerate) and 1 (normal).
constexpr double kOhcImpairment = 1.0;  // Normal.

// Inner Hair Cell (IHC) impairment coefficient. Should be between 0
// (degenerate) and 1 (normal).
constexpr double kIhcImpairment = 1.0;  // Normal.

// Minimum characteristic frequency in Hz (human).
constexpr double kCharFrequencyMin = 25.0;

// Maximum characteristic frequency in Hz (human).
constexpr double kCharFrequencyMax = 20E+3;

}  // namespace

void Zilany2014IhcModel::Init(const StimulusConfig &stimulus_config) {
  AuditoryModel::Init(stimulus_config);

  int num_channels = stimulus_config.num_channels();
  if (num_channels <= 0) {   // Use defaults.
    num_channels = kDefaultNumChannels;
  }
  mutable_stimulus_config()->set_num_channels(num_channels);

  // Compute linearly spaced characteristic frequencies.
  auto &cf_channels = channel_properties_.center_frequencies;
  cf_channels = erb::GetFrequenciesHz(
      num_channels, kCharFrequencyMin, kCharFrequencyMax);
}

ChannelProperties Zilany2014IhcModel::GetChannelProperties() const {
  return channel_properties_;
}

void Zilany2014IhcModel::ProcessSegment(
    const Eigen::Ref<const Eigen::ArrayXXd> &input, Eigen::ArrayXXd *output) {
  // Estimate the IHC transmembrane voltages (aka membrane potentials).
  const auto &cf_channels = channel_properties_.center_frequencies;
  const int sample_rate = stimulus_config().sample_rate();
  GOOGLE_CHECK_LE(100E3, sample_rate)
      << "Sample rate lower than lower bound: " << sample_rate;
  GOOGLE_CHECK_GE(500E3, sample_rate)
      << "Sample rate higher than upper bound: " << sample_rate;
  const Eigen::ArrayXd &signal = input.row(0);
  const int64_t num_samples = input.cols();
  const int num_channels = stimulus_config().num_channels();
  *output = Eigen::ArrayXXd::Zero(num_channels, num_samples);
  Eigen::RowArrayXd ihc_voltages(num_samples);
  for (int i = 0; i < num_channels; ++i) {
    IhcModel(signal,
             cf_channels[i],  // characteristic frequency
             1,               // number of repetitions
             1.0 / sample_rate, signal.size(), kOhcImpairment,
             kIhcImpairment, SPECIES_HUMAN_SHERA, &ihc_voltages);
    output->row(i).swap(ihc_voltages);
  }
}

}  // namespace audition
}  // namespace eidos

// Local Variables:
// mode: c++
// End:
