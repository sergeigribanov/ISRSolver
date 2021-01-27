#include "RangeInterpolator.hpp"

RangeInterpolator::RangeInterpolator(const std::tuple<bool, int, int>& rangeSettings,
                                     const Eigen::VectorXd& extCMEnergies) {
  // Vector extCMEnergies contains center-of-mass energies including threshold energy
  bool interpType;
  int rangeIndexMin;
  int rangeIndexMax;
  std::tie(interpType, rangeIndexMin, rangeIndexMax) = rangeSettings;
  _minEnergy = extCMEnergies(rangeIndexMin);
  _maxEnergy = extCMEnergies(rangeIndexMax + 1);
  if (rangeIndexMin > 0) {
    _beginIndex = rangeIndexMin;
  } else {
    _beginIndex = rangeIndexMin + 1;
  }
  const int numOfElements = rangeIndexMin > 0?
                            rangeIndexMax - rangeIndexMin + 2:
                            rangeIndexMax - rangeIndexMin + 1;
  _acc = std::vector<gsl_interp_accel*>(numOfElements);
  _spline = std::vector<gsl_spline*>(numOfElements);
  Eigen::VectorXd yi = Eigen::VectorXd::Zero(extCMEnergies.rows());
  for (int i = 0; i < numOfElements; ++i) {
    _acc[i] = gsl_interp_accel_alloc();
    if (interpType) {
      _spline[i] = gsl_spline_alloc(gsl_interp_cspline, extCMEnergies.rows());
    } else {
      _spline[i] = gsl_spline_alloc(gsl_interp_linear, extCMEnergies.rows());
    }
    yi(_beginIndex + i) = 1.;
    gsl_spline_init(_spline[i], extCMEnergies.data(), yi.data(), extCMEnergies.rows());
    yi(_beginIndex + i) = 0.;
  }
  _beginIndex--;
}

RangeInterpolator::~RangeInterpolator() {
  for (std::size_t i = 0; i < _spline.size(); ++i) {
    gsl_spline_free(_spline[i]);
    gsl_interp_accel_free(_acc[i]);
  }
}

double RangeInterpolator::basisEval(int csIndex, double energy) const {
  const int index = csIndex - _beginIndex;
  return gsl_spline_eval(_spline[index], energy, _acc[index]);
}

double RangeInterpolator::basisDerivEval(int csIndex, double energy) const {
  const int index = csIndex - _beginIndex;
  return gsl_spline_eval_deriv(_spline[index], energy, _acc[index]);
}

double RangeInterpolator::eval(const Eigen::VectorXd& y, double energy) const {
  double result = 0;
  int index;
  for (std::size_t i = 0; i < _spline.size(); ++i) {
    index = _beginIndex + i;
    result += basisEval(index, energy) * y(index);
  }
  return result;
}

double RangeInterpolator::derivEval(const Eigen::VectorXd& y, double energy) const {
  double result = 0;
  int index;
  for (std::size_t i = 0; i < _spline.size(); ++i) {
    index = _beginIndex + i;
    result += basisDerivEval(index, energy) * y(index);
  }
  return result;
}
