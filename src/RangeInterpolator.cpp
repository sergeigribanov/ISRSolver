#include <iostream>
#include "kuraev_fadin.hpp"
#include "RangeInterpolator.hpp"

RangeInterpolator::RangeInterpolator(const RangeInterpolator& rinterp):
    _cspline(rinterp._cspline),
    _beginIndex(rinterp._beginIndex), _minEnergy(rinterp._minEnergy),
    _maxEnergy(rinterp._maxEnergy), _acc(rinterp._acc),
    _spline(rinterp._spline) {}

RangeInterpolator::RangeInterpolator(const std::tuple<bool, int, int>& rangeSettings,
                                     const Eigen::VectorXd& extCMEnergies) {
  // Vector extCMEnergies contains center-of-mass energies including threshold energy
  int rangeIndexMin;
  int rangeIndexMax;
  std::tie(_cspline, rangeIndexMin, rangeIndexMax) = rangeSettings;
  if (rangeIndexMin > 0) {
    _beginIndex = rangeIndexMin - 1;
  } else {
    _beginIndex = 0;
  }
  _minEnergy = extCMEnergies(rangeIndexMin);
  _maxEnergy = extCMEnergies(rangeIndexMax + 1);
  const int numOfElements = rangeIndexMin > 0?
                            rangeIndexMax - rangeIndexMin + 2:
                            rangeIndexMax - rangeIndexMin + 1;
  _acc = std::vector<std::shared_ptr<gsl_interp_accel>>(numOfElements);
  _spline = std::vector<std::shared_ptr<gsl_spline>>(numOfElements);
  Eigen::VectorXd yi = Eigen::VectorXd::Zero(extCMEnergies.rows());
  for (int i = 0; i < numOfElements; ++i) {
    _acc[i] = std::shared_ptr<gsl_interp_accel>(
        gsl_interp_accel_alloc(),
        [](gsl_interp_accel* tacc)
        {gsl_interp_accel_free(tacc);});
    if (_cspline) {
      _spline[i] = std::shared_ptr<gsl_spline>(
          gsl_spline_alloc(gsl_interp_cspline, extCMEnergies.rows()),
          [](gsl_spline* tspline)
          {gsl_spline_free(tspline);});
    } else {
      _spline[i] = std::shared_ptr<gsl_spline>(
          gsl_spline_alloc(gsl_interp_linear, extCMEnergies.rows()),
          [](gsl_spline* tspline)
          {gsl_spline_free(tspline);});
    }
    yi(_beginIndex + i + 1) = 1.;
    gsl_spline_init(_spline[i].get(), extCMEnergies.data(), yi.data(), extCMEnergies.rows());
    yi(_beginIndex + i + 1) = 0.;
  }
}

RangeInterpolator::~RangeInterpolator() {}

double RangeInterpolator::basisEval(int csIndex, double energy) const {
  const int index = csIndex - _beginIndex;
  return gsl_spline_eval(_spline[index].get(), energy, _acc[index].get());
}

double RangeInterpolator::basisDerivEval(int csIndex, double energy) const {
  const int index = csIndex - _beginIndex;
  double koeff = 1.;
  if (!_cspline) {
    koeff = 1. - 1.e-12;
  }
  return gsl_spline_eval_deriv(_spline[index].get(), energy * koeff, _acc[index].get());
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

int RangeInterpolator::getBeginIndex() const {
  return _beginIndex;
}

double RangeInterpolator::getMinEnergy() const {
  return _minEnergy;
}

double RangeInterpolator::getMaxEnergy() const {
  return _maxEnergy;
}

bool RangeInterpolator::hasCSIndex(int csIndex) const {
  const int index = csIndex - _beginIndex;
  const int n = _spline.size();
  return (index >= 0) && (index < n);
}

bool RangeInterpolator::isEnergyInRange(double energy) const {
  return (energy > _minEnergy) && (energy <= _maxEnergy);
}

double RangeInterpolator::evalKuraevFadinBasisIntegral(
    int energyIndex, int csIndex, const std::function<double(double, double)>& efficiency) const {
  if (_cspline) {
    return _evalKuraevFadinIntegralCSpline(energyIndex, csIndex, efficiency);
  }
  return _evalKuraevFadinIntegralLinear(energyIndex, csIndex, efficiency);
}

double RangeInterpolator::_evalKuraevFadinIntegralLinear(
    int energyIndex, int csIndex, const std::function<double(double, double)>& efficiency) const {
  if (csIndex > energyIndex) {
    return 0;
  }
  const int index = csIndex - _beginIndex;
  std::function<double(double)> fcn0 =
      [](double energy) {
        return 1.;
      };
  std::function<double(double)> fcn1 =
      [](double energy) {
        return energy;
      };
  const double en = _spline[index].get()->x[energyIndex + 1];
  const double en0 = _spline[index].get()->x[csIndex];
  const double en1 = _spline[index].get()->x[csIndex + 1];
  const double x0 = 1 - std::pow(en0 / en, 2);
  const double x1 = 1 - std::pow(en1 / en, 2);
  const double i00 = kuraev_fadin_convolution(en, fcn0, x1, x0, efficiency);
  const double i01 = kuraev_fadin_convolution(en, fcn1, x1, x0, efficiency);
  double result = 1. / (en1 - en0) * (i01 - en0 * i00);
  if (csIndex + 2 < (int) _spline[index].get()->size) {
    const double en2 = _spline[index].get()->x[csIndex + 2];
    const double x2 = 1 - std::pow(std::min(en2, en) / en, 2);
    const double i10 = kuraev_fadin_convolution(en, fcn0, x2, x1, efficiency);
    const double i11 = kuraev_fadin_convolution(en, fcn1, x2, x1, efficiency);
    result += 1. / (en1 - en2) * (i11 - en1 * i10) + i10;
  }
  return result;
}

double RangeInterpolator::_evalKuraevFadinIntegralCSpline(
    int energyIndex, int csIndex, const std::function<double(double, double)>& efficiency) const {
  const int index = csIndex - _beginIndex;
  std::function<double(double)> fcn =
      [index, this] (double energy) {
        double result =  gsl_spline_eval(this->_spline[index].get(), energy, this->_acc[index].get());
        return result;
      };
  const double x_min = 1 - std::pow(_spline[index].get()->x[0] / _minEnergy, 2);
  const double x_max = 1 - std::pow(_spline[index].get()->x[0] / _spline[index].get()->x[energyIndex + 1], 2);
  double integralResult = kuraev_fadin_convolution(_spline[index].get()->x[energyIndex + 1], fcn, x_min, x_max, efficiency);
  return integralResult;
}
