#include "BaseRangeInterpolator.hpp"

BaseRangeInterpolator::BaseRangeInterpolator(
    int rangeIndexMin, int rangeIndexMax,
    const Eigen::VectorXd& extCMEnergies):
    _numberOfSegments(_evalNumOfSegments(rangeIndexMin, rangeIndexMax)),
    _beginIndex(_evalBeginIndex(rangeIndexMin, rangeIndexMax)),
    _minEnergy(extCMEnergies(rangeIndexMin)),
    _maxEnergy(extCMEnergies(rangeIndexMax + 1)) {}

BaseRangeInterpolator::BaseRangeInterpolator(const BaseRangeInterpolator& rinterp):
    _numberOfSegments(rinterp._numberOfSegments),
    _beginIndex(rinterp._beginIndex),
    _minEnergy(rinterp._minEnergy),
    _maxEnergy(rinterp._maxEnergy) {}

bool BaseRangeInterpolator::hasCSIndex(int csIndex) const {
  const int index = csIndex - _beginIndex;
  return (index >= 0) && (index < _numberOfSegments);
}

bool BaseRangeInterpolator::isEnergyInRange(double energy) const {
  return (energy > _minEnergy) && (energy <= _maxEnergy);
}

int BaseRangeInterpolator::getBeginIndex() const {
  return _beginIndex;
}

int BaseRangeInterpolator::getNumberOfSegments() const {
  return _numberOfSegments;
}

int BaseRangeInterpolator::_evalNumOfSegments(int rangeIndexMin,
                                              int rangeIndexMax) {
  const int numOfSegments = rangeIndexMin > 0?
                            rangeIndexMax - rangeIndexMin + 2:
                            rangeIndexMax - rangeIndexMin + 1;
  return numOfSegments;
}

int BaseRangeInterpolator::_evalBeginIndex(int rangeIndexMin,
                                           int rangeIndexMax) {
  const int beginIndex = rangeIndexMin > 0?
                         rangeIndexMin - 1:
                         0;
  return beginIndex;
}

double BaseRangeInterpolator::eval(const Eigen::VectorXd& y, double energy) const {
  double result = 0;
  for (int i = 0; i < _numberOfSegments; ++i) {
    const int index = _beginIndex + i;
    result += basisEval(index, energy) * y(index);
  }
  return result;
}

double BaseRangeInterpolator::derivEval(const Eigen::VectorXd& y, double energy) const {
  double result = 0;
  for (int i = 0; i < _numberOfSegments; ++i) {
    const int index = _beginIndex + i;
    result += basisDerivEval(index, energy) * y(index);
  }
  return result;
}

double BaseRangeInterpolator::getMinEnergy() const {
  return _minEnergy;
}

double BaseRangeInterpolator::getMaxEnergy() const {
  return _maxEnergy;
}
