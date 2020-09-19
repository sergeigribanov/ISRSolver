#define _USE_MATH_DEFINES
#include "kuraev_fadin.hpp"

#include <cmath>
#include <functional>

#include "integration.hpp"
#include "physical_constants.hpp"

double fLog(double s) { return std::log(s / ELECTRON_M / ELECTRON_M); }

double fBeta(double s) { return (2 * ALPHA_QED) / M_PI * (fLog(s) - 1); }

double kuraev_fadin_kernel(double x, double s) {
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

double kuraev_fadin_kernel_multiplication(
    double x, double s, const std::function<double(double)>& fcn) {
  return fcn(s * (1 - x)) * kuraev_fadin_kernel(x, s);
}

double kuraev_fadin_convolution(double s,
                                const std::function<double(double)>& fcn,
                                double min_x, double max_x) {
  std::function<double(double)> fcnConv = [s, &fcn](double x) {
    return kuraev_fadin_kernel_multiplication(x, s, fcn);
  };
  double error;
  double E = std::sqrt(s);
  double x0 = 2 * ELECTRON_M / E;
  double result;
  if (min_x < x0) {
    result = integrateS(fcnConv, min_x, x0, error) +
             integrate(fcnConv, x0, max_x, error);
  } else {
    result = integrate(fcnConv, min_x, max_x, error);
  }
  return result;
}
