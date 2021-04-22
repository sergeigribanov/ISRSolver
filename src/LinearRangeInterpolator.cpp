#include "Integration.hpp"
#include "KuraevFadin.hpp"
#include "LinearRangeInterpolator.hpp"

std::function<double(double)> LinearRangeInterpolator::_fcn0 =
    [](double) {return 1.;};
std::function<double(double)> LinearRangeInterpolator::_fcn1 =
    [](double energy) {return energy;};

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
  if (energy <= _extCMEnergies(csIndex)) {
    return 0;
  }
  if (energy <= _extCMEnergies(csIndex + 1)) {
    return _c01(csIndex);
  }
  if (csIndex + 2 == _extCMEnergies.rows()) {
    return 0;
  }
  if (energy <= _extCMEnergies(csIndex + 2)) {
    return _c11(csIndex);
  }
  return 0;
}

double LinearRangeInterpolator::evalKuraevFadinBasisIntegral(
    int energyIndex, int csIndex,
    const std::function<double(double, double)>& efficiency) const {
  if (csIndex > energyIndex) {
    return 0;
  }
  double result = _evalKuraevFadinBasisIntegralFirstTriangle(
      energyIndex, csIndex, efficiency);
  result += _evalKuraevFadinBasisIntegralSecondTriangle(
      energyIndex, csIndex, efficiency);
  return result;
}

double LinearRangeInterpolator::_evalKuraevFadinBasisIntegralFirstTriangle(
    int energyIndex, int csIndex,
    const std::function<double(double, double)>& efficiency) const {
  const double enc = _extCMEnergies(csIndex + 1);
  if (enc > _maxEnergy || enc <= _minEnergy) {
    return 0.;
  }
  const double en = _extCMEnergies(energyIndex + 1);
  const double x0 = std::max(0., 1 - std::pow(enc / en, 2));
  const double x1 = 1 - std::pow(_extCMEnergies(csIndex) / en, 2);
  const double i00 = kuraev_fadin_convolution(en, _fcn0, x0, x1, efficiency);
  const double i01 = kuraev_fadin_convolution(en, _fcn1, x0, x1, efficiency);
  return _c01(csIndex) * i01 + _c00(csIndex) * i00;
}

double LinearRangeInterpolator::_evalKuraevFadinBasisIntegralSecondTriangle(
    int energyIndex, int csIndex,
    const std::function<double(double, double)>& efficiency) const {
  if (csIndex + 2 >= _extCMEnergies.rows()) {
    return 0.;
  }
  const double encp = _extCMEnergies(csIndex + 2);
  if (encp > _maxEnergy || encp <= _minEnergy) {
    return 0.;
  }
  const double en = _extCMEnergies(energyIndex + 1);
  const double x0 = std::max(0., 1 - std::pow(encp / en, 2));
  const double x1 = 1 - std::pow(_extCMEnergies(csIndex + 1) / en, 2);
  const double i10 = kuraev_fadin_convolution(en, _fcn0, x0, x1, efficiency);
  const double i11 = kuraev_fadin_convolution(en, _fcn1, x0, x1, efficiency);
  return _c11(csIndex) * i11 + _c10(csIndex) * i10;
}

double LinearRangeInterpolator::_evalIntegralBasisFirstTriangle(int csIndex) const {
  const double enc = _extCMEnergies(csIndex + 1);
  if (enc > _maxEnergy || enc <= _minEnergy) {
    return 0.;
  }
  double error;
  const double i00 = integrate(_fcn0, _extCMEnergies(csIndex), enc, error);
  const double i01 = integrate(_fcn1, _extCMEnergies(csIndex), enc, error);
  return _c01(csIndex) * i01 + _c00(csIndex) * i00;
}

double LinearRangeInterpolator::_evalIntegralBasisSecondTriangle(int csIndex) const {
  if (csIndex + 2 >= _extCMEnergies.rows()) {
    return 0.;
  }
  const double encp = _extCMEnergies(csIndex + 2);
  if (encp > _maxEnergy || encp <= _minEnergy) {
    return 0.;
  }
  double error;
  const double enc = _extCMEnergies(csIndex + 1);
  const double i10 = integrate(_fcn0, enc, encp, error);
  const double i11 = integrate(_fcn1, enc, encp, error);
  return _c11(csIndex) * i11 + _c10(csIndex) * i10;
}

double LinearRangeInterpolator::evalIntegralBasis(int csIndex) const {
  double result = _evalIntegralBasisFirstTriangle(csIndex);
  result += _evalIntegralBasisSecondTriangle(csIndex);
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
