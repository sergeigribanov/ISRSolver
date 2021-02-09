#include "kuraev_fadin.hpp"
#include "LinearRangeInterpolator.hpp"

#include <iostream>

// Vector extCMEnergies contains center-of-mass energies including threshold energy
// TO DO: replace _extCMEnergies with pointer
LinearRangeInterpolator::LinearRangeInterpolator(
    int rangeIndexMin, int rangeIndexMax,
    const Eigen::VectorXd& extCMEnergies):
    BaseRangeInterpolator(rangeIndexMin, rangeIndexMax, extCMEnergies),
    _extCMEnergies(extCMEnergies) {}

LinearRangeInterpolator::~LinearRangeInterpolator() {}

LinearRangeInterpolator::LinearRangeInterpolator(const LinearRangeInterpolator& rinterp):
    BaseRangeInterpolator(rinterp),
    _extCMEnergies(rinterp._extCMEnergies) {}

double LinearRangeInterpolator::basisEval(int csIndex, double energy) const {
  if (energy <= _extCMEnergies(csIndex)) {
    return 0;
  }
  if (energy <= _extCMEnergies(csIndex + 1)) {
    return _c01(csIndex) * energy + _c00(csIndex);
  }
  if (csIndex + 2 == _extCMEnergies.rows()) {
    return 1;
  }
  if (energy <= _extCMEnergies(csIndex + 2)) {
    return _c11(csIndex) * energy + _c10(csIndex);
  }
  return 0;
}

double LinearRangeInterpolator::basisDerivEval(int csIndex, double energy) const {
  if (energy <= _extCMEnergies(csIndex) ||
      csIndex + 2 == _extCMEnergies.rows()) {
    return 0;
  }
  if (energy <= _extCMEnergies(csIndex + 1)) {
    return _c01(csIndex);
  }
  return _c11(csIndex);
}

double LinearRangeInterpolator::evalKuraevFadinBasisIntegral(
    int energyIndex, int csIndex,
    const std::function<double(double, double)>& efficiency) const {
  if (csIndex > energyIndex) {
    return 0;
  }
  std::function<double(double)> fcn0 =
      [](double energy) {
        return 1.;
      };
  std::function<double(double)> fcn1 =
      [](double energy) {
        return energy;
      };
  const double en = _extCMEnergies(energyIndex + 1);
  const double x1 = 1 - std::pow(_extCMEnergies(csIndex + 1) / en, 2);
  const double x2 = 1 - std::pow(_extCMEnergies(csIndex) / en, 2);
  const double i00 = kuraev_fadin_convolution(en, fcn0, x1, x2, efficiency);
  const double i01 = kuraev_fadin_convolution(en, fcn1, x1, x2, efficiency);
  double result = _c01(csIndex) * i01 + _c00(csIndex) * i00;
  if (csIndex + 2 < _extCMEnergies.rows()) {
    const double x0 = std::max(0., 1 - std::pow(_extCMEnergies(csIndex + 2) / en, 2));
    const double i10 = kuraev_fadin_convolution(en, fcn0, x0, x1, efficiency);
    const double i11 = kuraev_fadin_convolution(en, fcn1, x0, x1, efficiency);
    result += _c11(csIndex) * i11 + _c10(csIndex) * i10;
  }
  return result;
}

double LinearRangeInterpolator::_c00(int csIndex) const {
  return _extCMEnergies(csIndex) /
      (_extCMEnergies(csIndex) - _extCMEnergies(csIndex + 1));
}

double LinearRangeInterpolator::_c01(int csIndex) const {
  return 1. / (_extCMEnergies(csIndex + 1) - _extCMEnergies(csIndex));
}

double LinearRangeInterpolator::_c10(int csIndex) const {
  return 1. + _extCMEnergies(csIndex + 1) / (_extCMEnergies(csIndex + 2) - _extCMEnergies(csIndex + 1));
}

double LinearRangeInterpolator::_c11(int csIndex) const {
  return 1. / (_extCMEnergies(csIndex + 1) - _extCMEnergies(csIndex + 2));
}
