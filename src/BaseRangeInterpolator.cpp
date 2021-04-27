#include "BaseRangeInterpolator.hpp"

/**
 * Constructor
 * @param rangeIndexMin a minimum index of sub range (a range between
 two neighboring center-of-mass energy points) that used in interpolation
 range
 * @param rangeIndexMax a maximum index of sub range that used in
 interpolation range
 * @param extCMEnergies extended vector of center-of-mass energies (contains
 threshold energy)
*/
BaseRangeInterpolator::BaseRangeInterpolator(
    int rangeIndexMin, int rangeIndexMax,
    const Eigen::VectorXd& extCMEnergies):
    _numberOfSegments(_evalNumOfSegments(rangeIndexMin, rangeIndexMax)),
    _beginIndex(_evalBeginIndex(rangeIndexMin, rangeIndexMax)),
    _minEnergy(extCMEnergies(rangeIndexMin)),
    _maxEnergy(extCMEnergies(rangeIndexMax + 1)) {}

/**
 * Copy constructor
 */
BaseRangeInterpolator::BaseRangeInterpolator(const BaseRangeInterpolator& rinterp):
    _numberOfSegments(rinterp._numberOfSegments),
    _beginIndex(rinterp._beginIndex),
    _minEnergy(rinterp._minEnergy),
    _maxEnergy(rinterp._maxEnergy) {}

/**
 * This method returns true if cross section point
 * with csIndex belongs interpolation sub range
 * @param csIndex an index of a cross section point
 */
bool BaseRangeInterpolator::hasCSIndex(int csIndex) const {
  const int index = csIndex - _beginIndex;
  return (index >= 0) && (index < _numberOfSegments);
}

/**
 * This method returns true if energy belongs
 to the considered interpolation sub range
 * @param energy a center-of-mass energy
 */
bool BaseRangeInterpolator::isEnergyInRange(double energy) const {
  return (energy > _minEnergy) && (energy <= _maxEnergy);
}

/**
 * This method returns the index of the first segment
 (a segment with lowest center-of-mass energies) that belongs
 to the considered sub range
*/
int BaseRangeInterpolator::getBeginIndex() const {
  return _beginIndex;
}

/**
 * This method returns a number of center-of-mass energy
 segments (from one
 cross section point to another) in the considered
 sub range.
 * The segment indices are defined as follows.
 * If index=0, the center-of-mass energy segment
 between threshold energy and the first center-of-mass
 energy is considered. In the case when index=1,
 the center-of-mass energy segment between the first
 and the second center-of-mass energy points is considered.
*/
int BaseRangeInterpolator::getNumberOfSegments() const {
  return _numberOfSegments;
}

/**
 * Evaluating number of segments of the current
 interpolation sub range
 * @param rangeIndexMin a minimum sub range index
 * @param rangeIndexMax a maximum sub range index
 */
int BaseRangeInterpolator::_evalNumOfSegments(
    int rangeIndexMin,
    int rangeIndexMax) {
  const int numOfSegments = rangeIndexMin > 0?
                            rangeIndexMax - rangeIndexMin + 2:
                            rangeIndexMax - rangeIndexMin + 1;
  return numOfSegments;
}

/**
 * Evaluating the begin index
 * @param rangeIndexMin a minimum sub range index
 * @param rangeIndexMax a maximum sub range index
 */
int BaseRangeInterpolator::_evalBeginIndex(int rangeIndexMin,
                                           int rangeIndexMax) {
  const int beginIndex = rangeIndexMin > 0?
                         rangeIndexMin - 1:
                         0;
  return beginIndex;
}

/**
 * This method is used to eval interpolation
 value at the certain center-of-mass energy
 * @param y a data
 * @param energy a center-of-mass energy
 */
double BaseRangeInterpolator::eval(const Eigen::VectorXd& y, double energy) const {
  double result = 0;
  for (int i = 0; i < _numberOfSegments; ++i) {
    const int index = _beginIndex + i;
    /**
     * Adding the contribution of the i-th center-of-mass energy point
     */
    result += basisEval(index, energy) * y(index);
  }
  return result;
}

/**
 * This method is used to eval interpolation
 derivative value at the certain center-of-mass energy
 * @param y a data
 * @param energy a center-of-mass energy
 */
double BaseRangeInterpolator::derivEval(const Eigen::VectorXd& y, double energy) const {
  double result = 0;
  for (int i = 0; i < _numberOfSegments; ++i) {
    const int index = _beginIndex + i;
    /**
     * Adding the contribution of the i-th center-of-mass energy point
    */
    result += basisDerivEval(index, energy) * y(index);
  }
  return result;
}

/**
 * This method is used to get minimum energy
 */
double BaseRangeInterpolator::getMinEnergy() const {
  return _minEnergy;
}

/**
 * This method is used to get maximum energy
 */
double BaseRangeInterpolator::getMaxEnergy() const {
  return _maxEnergy;
}
