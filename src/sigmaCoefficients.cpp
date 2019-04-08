#include "sigmaCoefficients.h"
#include "radiatorFadinKuraev.h"

double getDet2 (double a, double b,
		double c, double d) {
  return a * d - b * c;
}

std::vector<double> getLinearSigmaCoeffs(double xm, double xi,
					 double s, double min_x, double max_x) {
  double det2 = getDet2 (xm, 1, xi, 1);
  double integral0 = radIntegral (s, min_x, max_x, 0);
  double integral1 = radIntegral (s, min_x, max_x, 1);
  double cm = (integral1 - xi * integral0) / det2;
  double ci = (-integral1 + xm * integral0) / det2;
  return std::vector<double> {cm, ci};
}
