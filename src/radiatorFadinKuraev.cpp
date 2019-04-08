#include "radiatorFadinKuraev.h"
#include "phyConstants.h"
#include "getIntegral.h"
#include <math.h>
#include <functional>

double fLog (double s) {
  return log(s / ELECTRON_M / ELECTRON_M);
}

double fBeta (double s) {
  return (2 * ALPHA_QED) / M_PI * (fLog (s) - 1);
}

double radiatorFadinKuraev (double x, double s) {
  double lnX = log(x);
  double mX = 1 - x;
  double lnMX = log(mX);
  double sM = s / ELECTRON_M / ELECTRON_M;
  double logSM = log(sM);
  double E = 0.5 * sqrt(s);
  double res = 0;
    
  double part1 = fBeta (s) * pow(x, fBeta(s) - 1) *
    (1 + ALPHA_QED / M_PI * (M_PI * M_PI / 3 - 0.5) +
     0.75 * fBeta (s) - fBeta (s) * fBeta (s) / 24 *
     (fLog (s) / 3 + 2 * M_PI * M_PI - 9.25));

  double part2 = -fBeta (s) * (1 - 0.5 * x);

  double part3 = 0.125 * fBeta (s) * fBeta (s) *
    (-4 * (2 - x) * lnX - (1 + 3 * mX * mX) / x * lnMX -
     6 + x);

  res = part1 + part2 + part3;

  if (x > 2 * ELECTRON_M / E) {
    
    double subsubpart1 = 2 * lnX + logSM - 5. / 3;
  
    double subpart1 = pow(x - 2 * ELECTRON_M / E , fBeta (s) ) / 6 / x *
      pow(subsubpart1, 2) * (2 - 2 * x + x * x + fBeta (s) / 3 * subsubpart1);
    
    double subpart2 = 0.5 * fLog (s) * fLog (s) * (2. / 3 * (1 - mX * mX *mX) / mX +
						 (2 - x) * lnMX + 0.5 * x);
    res += (ALPHA_QED / M_PI) * (ALPHA_QED / M_PI) * (subpart1 + subpart2);
  }

  return res;
}

double polRadiator (double x, double s, int n) {
  if (n == 0) {
    return radiatorFadinKuraev(x, s);
  } else if (n == 1) {
    return x * radiatorFadinKuraev(x, s);
  } else if (n == 2) {
    return x * x * radiatorFadinKuraev(x, s);
  } else {
    return pow(x, n) * radiatorFadinKuraev(x, s);
  }
}

double radIntegral (double s, double min_x, double max_x, int n) {
  std::function<double (double) > fcn = [s, n] (double x) {
    return polRadiator (x, s, n);
  };
  double error;
  double part1;
  double part2;
  double E = sqrt(s);
  double x1 = 2 * ELECTRON_M / E;
  if (min_x < x1) {
    part1 = getIntegralS(fcn, min_x, x1, error);
    part2 = getIntegral(fcn, x1, max_x, error);
  } else {
    part1 = 0;
    part2 = getIntegral(fcn, min_x, max_x, error);
  }
  return part1 + part2;
}

double getFadinConv(double x,
                    double s,
                    const std::function<double(double)>& fcn) {
  return fcn(s * (1 - x)) * radiatorFadinKuraev(x, s);
}

double getFadinIntegral(double s,
                        const std::function<double(double)>& fcn,
                        double max_x) {
  std::function<double(double)> fcnConv = [s, &fcn] (double x) {
    return getFadinConv(x, s, fcn);
  };
  double error;
  double E = sqrt(s);
  double part1 = getIntegralS(fcnConv, 0, 2 * ELECTRON_M / E, error);
  double part2 = getIntegral(fcnConv, 2 * ELECTRON_M / E, max_x, error);
  return part1 + part2;
}
