#define _USE_MATH_DEFINES

#include <cmath>
#include <functional>

#include "ISRSolver/Integration.hpp"
#include "ISRSolver/PhysicalConstants.hpp"
#include "ISRSolver/KuraevFadin.hpp"

double fLog(double s) { return std::log(s / ELECTRON_M / ELECTRON_M); }

double fBeta(double s) { return (2 * ALPHA_QED) / M_PI * (fLog(s) - 1); }

double kernelKuraevFadin(double x, double s) {
  double lnX = std::log(x);
  double mX = 1 - x;
  double lnMX = std::log(mX);
  double sM = s / ELECTRON_M / ELECTRON_M;
  double logSM = std::log(sM);
  double E = 0.5 * std::sqrt(s);
  double res = 0;

  double part1 =
      fBeta(s) * std::pow(x, fBeta(s) - 1) *
      (1 + ALPHA_QED / M_PI * (M_PI * M_PI / 3 - 0.5) + 0.75 * fBeta(s) -
       fBeta(s) * fBeta(s) / 24 * (fLog(s) / 3 + 2 * M_PI * M_PI - 9.25));

  double part2 = -fBeta(s) * (1 - 0.5 * x);

  double part3 = 0.125 * fBeta(s) * fBeta(s) *
                 (-4 * (2 - x) * lnX - (1 + 3 * mX * mX) / x * lnMX - 6 + x);

  res = part1 + part2 + part3;

  if (x > 2 * ELECTRON_M / E) {
    double subsubpart1 = 2 * lnX + logSM - 5. / 3;

    double subpart1 = std::pow(x - 2 * ELECTRON_M / E, fBeta(s)) / 6 / x *
                      std::pow(subsubpart1, 2) *
                      (2 - 2 * x + x * x + fBeta(s) / 3 * subsubpart1);

    double subpart2 =
        0.5 * fLog(s) * fLog(s) *
        (2. / 3 * (1 - mX * mX * mX) / mX + (2 - x) * lnMX + 0.5 * x);
    res += (ALPHA_QED / M_PI) * (ALPHA_QED / M_PI) * (subpart1 + subpart2);
  }

  return res;
}

double kernelMultiplicationKuraevFadin(
    double x, double energy, const std::function<double(double)>& fcn) {
  return fcn(energy * std::sqrt(1 - x)) * kernelKuraevFadin(x, energy * energy);
}

double convolutionKuraevFadin(double energy,
                                const std::function<double(double)>& fcn,
                                double min_x, double max_x,
				const std::function<double(double, double)>& efficiency) {
  std::function<double(double)> fcnConv = [energy, &fcn, &efficiency](double x) {
    return kernelMultiplicationKuraevFadin(x, energy, fcn) * efficiency(x, energy);
  };
  double error;
  double x0 = 4 * ELECTRON_M / energy;
  double result;
  if (min_x < x0) {
    result = integrateS(fcnConv, min_x, x0, error) +
             integrate(fcnConv, x0, max_x, error);
  } else {
    result = integrate(fcnConv, min_x, max_x, error);
  }
  return result;
}

double convolutionKuraevFadin_1d(double energy,
                                 const std::function<double(double)>& fcn,
                                 double min_x, double max_x,
                                 const std::function<double(double)>& efficiency) {
  std::function<double(double)> fcnConv = [energy, &fcn, &efficiency](double x) {
    return kernelMultiplicationKuraevFadin(x, energy, fcn) * efficiency(x);
  };
  double error;
  double x0 = 4 * ELECTRON_M / energy;
  double result;
  if (min_x < x0) {
    result = integrateS(fcnConv, min_x, x0, error) +
             integrate(fcnConv, x0, max_x, error);
  } else {
    result = integrate(fcnConv, min_x, max_x, error);
  }
  return result;
}


// mX = 1 - x
double funcDGamma(double mX, double s) {
  const double fL = fLog(s);
  const double beta = fBeta(s);
  const double x = 1. - mX;
  const double lnX = std::log(x);
  const double lnMX = std::log(mX);
  return 0.5 * beta * std::pow(mX, 0.5 * beta - 1.) *
      (1. + 3. / 8 * beta - 1. / 48 * beta * beta * (1. / 3 * fL + M_PI * M_PI - 47. / 8))
      - 0.25 * beta * (1. + x) - 1. / 32 * beta * beta *
      (4. * (1. + x) * lnMX + (1. + 3. * std::pow(x, 2.)) / mX * lnX + 5. + x);
}

// mX = 1 - x
double funcDEE(double mX, double s) {
  const double E = 0.5 * std::sqrt(s);
  if (mX - 2. * ELECTRON_M / E < 0) {
    return 0.;
  }
  const double fL = fLog(s);
  const double beta = fBeta(s);
  const double x = 1. - mX;
  const double sL = fL + 2. * std::log(mX);
  const double lnX = std::log(x);
  return std::pow(ALPHA_QED / M_PI, 2.) * (
      1. / 12 / mX * std::pow(mX - 2. * ELECTRON_M / E, 0.5 * beta) *
      std::pow(sL - 5. / 3, 2) * (1. + std::pow(x, 2.) + 1. / 6 * beta * (sL - 5. / 3))
      + 0.25 * fL * fL * (2. / 3 * (1. - std::pow(x, 3)) / x + 0.5 * mX + (1. + x) * lnX));
}

// mX = 1 - x
double funcD(double mX, double s) {
  return funcDGamma(mX, s) + funcDEE(mX, s);
}

double funcDIntegral(double mX, double s, double s_threshold) {
  const double mX_max = 1. - s_threshold / s;
  if (mX >= mX_max) {
    return 0.;
  }
  std::function<double(double)> fcn = [s, mX_max](double mY) {
    if (mY >= mX_max) {
      return 0.;
    }
    return funcD(mY, s);
  };
  double error;
  const double E = 0.5 * std::sqrt(s);
  const double mX0 = 2. * ELECTRON_M / E;
  if (mX < mX0) {
    return integrateS(fcn, mX, mX0, error) +
        integrate(fcn, mX0, mX_max, error);
  }
  return integrate(fcn, mX, mX_max, error);
}

double funcFIntegral(double x_min, double x_max, double s) {
  std::function<double(double)> fcn = [s](double y) {
    return kernelKuraevFadin(y, s);
  };
  double error;
  const double E = 0.5 * std::sqrt(s);
  const double x0 = 2. * ELECTRON_M / E;
  if (x_min < x0) {
    return integrateS(fcn, x_min, x0, error) +
        integrate(fcn, x0, x_max, error);
  }
  return integrate(fcn, x_min, x_max, error);
}

// mX = 1 - x
double majorantDFunction(double mX, double s, double koeff) {
  const double beta = fBeta(s);
  return koeff * 0.5 * beta * std::pow(mX, 0.5 * beta - 1.);
}
