#include <algorithm>
#include <fstream>
#include "Interpolator.hpp"
#include "LinearRangeInterpolator.hpp"
#include "CSplineRangeInterpolator.hpp"

/**
 * Default constructor
 */
Interpolator::Interpolator() {}

/**
 * Copy constructor
 */
Interpolator::Interpolator(const Interpolator& interp):
    _rangeInterpolators(interp._rangeInterpolators) {}

/**
 * Constructor
 * @param cmEnergies center-of-mass energies
 * @param thresholdEnergy a threshold energy
 */
Interpolator::Interpolator(const Eigen::VectorXd& cmEnergies,
                           double thresholdEnergy):
    Interpolator(defaultInterpRangeSettings(cmEnergies.rows()),
                 cmEnergies, thresholdEnergy) {}

/**
 * Constructor
 * @param interpRangeSettings an interpolation settings
 * @param cmEnergies center-of-mass energies
 * @param thresholdEnergy a threshold energy
 */
Interpolator::Interpolator(
    const std::vector<std::tuple<bool, int, int>>& interpRangeSettings,
    const Eigen::VectorXd& cmEnergies, double thresholdEnergy) {
  Eigen::VectorXd extCMEnergies(cmEnergies.rows() + 1);
  extCMEnergies.tail(cmEnergies.rows()) = cmEnergies;
  extCMEnergies(0) = thresholdEnergy;
  auto interpRSsorted = interpRangeSettings;
  /**
   * Sorting interpolation range settings
   */
  std::sort(interpRSsorted.begin(), interpRSsorted.end(),
            [](const std::tuple<bool, int, int>& t1,
               const std::tuple<bool, int, int>& t2) {
              return std::get<1>(t1) < std::get<1>(t2);});
  _rangeInterpolators.reserve(interpRSsorted.size());
  /**
   * Verifying interpolation settings
   */
  if (checkRangeInterpSettings(interpRSsorted,
                               cmEnergies)) {
    /**
     * Throw exception if interpolation settings are wrong
     */
    InterpRangeException ex;
    throw ex;
  }
  /**
   * Filling interpolators
   */
  for (const auto& el : interpRSsorted) {
    bool interpType;
    double rangeIndexMin;
    double rangeIndexMax;
    std::tie(interpType, rangeIndexMin, rangeIndexMax) = el;
    if (interpType) {
      /**
       * Filling cubic spline interpolators
       */
      _rangeInterpolators.push_back(
          std::shared_ptr<BaseRangeInterpolator>(new CSplineRangeInterpolator(
              rangeIndexMin, rangeIndexMax, extCMEnergies)));
    } else {
      /**
       * Filling piecewise linear interpolators
       */
      _rangeInterpolators.push_back(std::shared_ptr<BaseRangeInterpolator>(new LinearRangeInterpolator(
          rangeIndexMin, rangeIndexMax, extCMEnergies)));
    }
  }
}

/**
 * Check validity of range interpolation settings
 * @param sortedInterpRangeSettings an interpolation settings
 * @param cmEnergies a vector of center-of-mass energies
 */
bool Interpolator::checkRangeInterpSettings(
    const std::vector<std::tuple<bool, int, int>>&
    sortedInterpRangeSettings,
    const Eigen::VectorXd& cmEnergies) {
  /**
   * Returns true if interpolation settings are wrong
   */
  if(std::get<1>(*sortedInterpRangeSettings.begin()) != 0 ||
     std::get<2>(sortedInterpRangeSettings.back()) + 1
     != cmEnergies.rows()) {
    /**
     * (*sortedInterpRangeSettings.begin()) is the first interpolation range.
     * std::get<1>(*sortedInterpRangeSettings.begin()) is the index of the
     * first center-of-mass energy segment in the first interpolation range.
     * If this index is not equal to zero, interpolation settings are wrong.
     * std::get<2>(sortedInterpRangeSettings.back()) is the last segment of
     * the last interpolation ranges. The interpolation settings are wrong
     * if this index is not equal to the number of center-of-mass energy points
     * minus one.
     */
    return true;
  }
  /**
   * The last center-of-mass energy segment index in the previous
   * interpolation range
   */
  int pimax = -1;
  /**
   * Loop over interpolation ranges
   */
  for (const auto& el : sortedInterpRangeSettings) {
    /**
     * The first center-of-mass energy segment index in the current
     * interpolation range
     */
    int imin = std::get<1>(el);
    /**
     * The last center-of-mass energy segment index in the current
     * interpolation range
     */
    int imax = std::get<2>(el);
    if (imin > imax) {
      /**
       * Interpolation settings are wrong if imin > imax
       */
      return true;
    }
    if (pimax != -1 && pimax + 1 != imin) {
      /**
       * Interpolation settings are wrong if pimax + 1 = imin
       * (missed or taken twice center-of-mass energy segments)
       */
      return true;
    }
    pimax = imax;
  }
  return false;
}

/**
 * Constructor
 * @param pathToJSON a path to .json file with interpolation settings
 * @param cmEnergies a vector of center-of-mass energies
 * @param thresholdEnergy a threshold energy
 */
Interpolator::Interpolator(const std::string& pathToJSON,
                          const Eigen::VectorXd& cmEnergies,
                           double thresholdEnergy):
    Interpolator(Interpolator::loadInterpRangeSettings(pathToJSON),
                 cmEnergies, thresholdEnergy) {}

/**
 * Load interpolation settings from file
 * @param pathToJSON a path to a .json file with settings
 */
std::vector<std::tuple<bool, int, int>> Interpolator::loadInterpRangeSettings(
    const std::string& pathToJSON) {
  /**
   * Opening .json file with interpolation settings
   */
  std::ifstream fl(pathToJSON);
  json s;
  /**
   * Reading .json file to json object
   */
  fl >> s;
  std::vector<std::tuple<bool, int, int>> result;
  result.reserve(s.size());
  /**
   * Converting JSON data to appropriate format
   */
  for (const auto& el : s) {
    result.push_back(std::make_tuple(el[0].get<bool>(),
                                     el[1].get<int>(),
                                     el[2].get<int>()));
  }
  fl.close();
  return result;
}

/**
 * Destructor
 */
Interpolator::~Interpolator() {}

/**
 * Evaluating basis interpolation function
 */
double Interpolator::basisEval(int csIndex, double energy) const {
  if (energy <= getMinEnergy()) {
    /**
     * Returns zero if energy is less than minimum energy. The minimum
     * energy is a threshold energy in the case of interpolation.
     */
    return 0;
  }
  if (energy > getMaxEnergy()) {
    /**
     * If energy is greater than maximum energy returns interpolation
     * function value at point with the maximum center-of-mass energy
     */
    if (_rangeInterpolators.back().get()->hasCSIndex(csIndex)) {
      /**
       * If csIndex belongs to the last range interpolator, the interpolation
       * function of this interpolator is evaluated at the maximum
       * center-of-mass energy
       */
      return _rangeInterpolators.back().get()->basisEval(csIndex, getMaxEnergy());
    } else {
      /**
       * If csIndex doesn't belong to the last range interpolator, zero is
       * returned
       */
      return 0;
    }
  }
  double result = 0;
  /**
   * Loop over all interpolation ranges
   */
  for (const auto& rinterp : _rangeInterpolators) {
    if (rinterp.get()->hasCSIndex(csIndex) && rinterp.get()->isEnergyInRange(energy)) {
      /**
       * Adding contribution of the current interpolation range
       */
      result += rinterp.get()->basisEval(csIndex, energy);
    }
  }
  return result;
}

/**
 * Evaluating derivative of a basis interpolation function
 * @param csIndex an index of a corresponding cross section
 point
 @param energy is a center-of-mass energy
*/
double Interpolator::basisDerivEval(int csIndex, double energy) const {
  if (energy <= getMinEnergy() || energy > getMaxEnergy()) {
    /**
     * Derivative is zero if energy is out of limits
     */
    return 0;
  }
  double result = 0;
  /**
   * Loop over all interpolation ranges
   */
  for (const auto& rinterp : _rangeInterpolators) {
    if (rinterp.get()->hasCSIndex(csIndex) &&
        rinterp.get()->isEnergyInRange(energy)) {
      /**
       * Adding contribution of the current interpolation range
       */
      result += rinterp.get()->basisDerivEval(csIndex, energy);
    }
  }
  return result;
}

/**
 * Evaluating interpolation y data value at a certain
 center-of-mass energy
 * @param y a data
 * @param energy a center-of-mass energy
 */
double Interpolator::eval(const Eigen::VectorXd& y, double energy) const {
  if (energy <= getMinEnergy()) {
    /**
     * Returns zero if energy is less than minimum energy. The minimum
     * energy is a threshold energy in the case of interpolation.
     */
    return 0;
  }
  if (energy > getMaxEnergy()) {
    /**
     * If energy is greater than maximum energy returns interpolation
     * function value at point with the maximum center-of-mass energy
     */
    return _rangeInterpolators.back().get()->eval(y, getMaxEnergy());
  }
  double result = 0;
  /**
   * Loop over all interpolation ranges
   */
  for (const auto& rinterp : _rangeInterpolators) {
    if (rinterp.get()->isEnergyInRange(energy)) {
      /**
       * Adding contribution of the current interpolation range
       */
      result += rinterp.get()->eval(y, energy);
    }
  }
  return result;
}

/**
 * Evaluating interpolation derivative
 * @param y a data
 * @param energy a center-of-mass energy
 */
double Interpolator::derivEval(const Eigen::VectorXd& y, double energy) const {
  if (energy <= getMinEnergy() || energy > getMaxEnergy()) {
    /**
     * Derivative is zero if energy is out of limits
     */
    return 0;
  }
  double result = 0;
  /**
   * Loop over all interpolation ranges
   */
  for (const auto& rinterp : _rangeInterpolators) {
    if (rinterp.get()->isEnergyInRange(energy)) {
      /**
       * Adding contribution of the current interpolation range
       */
      result += rinterp.get()->derivEval(y, energy);
    }
  }
  return result;
}

/**
 * Use default interpolation settings (piecewise linear interpolation)
 * @param n a number of cross section points
 */
std::vector<std::tuple<bool, int, int>>
Interpolator::defaultInterpRangeSettings(int n) {
  return std::vector<std::tuple<bool, int, int>>(
      1, std::tuple<bool, int, int>(false, 0, n - 1));
}

/**
 * Getting the minimum center-of-mass energy associated with
 the considered interpolator
*/
double Interpolator::getMinEnergy() const {
  return _rangeInterpolators.begin()->get()->getMinEnergy();
}

/**
 * Getting the maximum center-of-mass energy associated with
 the considered interpolator
*/
double Interpolator::getMaxEnergy() const {
  return _rangeInterpolators.back().get()->getMaxEnergy();
}

/**
 * Getting the minimum center-of-mass energy
 that belongs interpolation range, which contains
 the center-of-mass energy segment with index csIndex
 * @param csIndex an index of center-of-mass energy segment
 */
double Interpolator::getMinEnergy(int csIndex) const {
  double result = _rangeInterpolators.back().get()->getMinEnergy();
  /**
   * Loop over all interpolation ranges in order to find
   * which interpolation range contains the center-of-mass
   * energy segment with the index csIndex
   */
  for (const auto& rinterp : _rangeInterpolators) {
    if(rinterp.get()->hasCSIndex(csIndex)) {
      /**
       * The corresponding interpolation range is found. Obtaining
       * the minimum center-of-mass energy.
       */
      result = std::min(result, rinterp.get()->getMinEnergy());
    }
  }
  return result;
}

/**
 * Getting the maximum center-of-mass energy
 * that belongs interpolation range, which contains
 * the center-of-mass energy segment with index csIndex
 * @param csIndex an index of center-of-mass energy segment
 */
double Interpolator::getMaxEnergy(int csIndex) const {
  double result = _rangeInterpolators[0].get()->getMaxEnergy();
  /**
   * Loop over all interpolation ranges in order to find
   * which interpolation range contains the center-of-mass
   * energy segment with the index csIndex
   */
  for (const auto& rinterp : _rangeInterpolators) {
    if(rinterp.get()->hasCSIndex(csIndex)) {
      /**
       * The corresponding interpolation range is found. Obtaining
       * the maximum center-of-mass energy.
       */
      result = std::max(result, rinterp.get()->getMaxEnergy());
    }
  }
  return result;
}

/**
 * Kuraev-Fadin convolution with a basis interpolation function
 * that corresponds to the cross section point indx csIndex.
 * @param energyIndex an index of center-of-mass energy at
 * which convolution is evaluated
 * @param csIndex an index of cross-section point that corresponds
 * to considered basis interpolation function
 * @param efficiency a detection efficiency
 */
double Interpolator::evalKuraevFadinBasisIntegral(
    int energyIndex, int csIndex,
    const std::function<double(double, double)>& efficiency) const {
  double result = 0;
  /**
   * Loop over all interpolation ranges
   */
  for (const auto& rinterp : _rangeInterpolators) {
    if(rinterp.get()->hasCSIndex(csIndex)) {
      /**
       * Evaluate contribution of each interpolation range
       */
      result += rinterp.get()->evalKuraevFadinBasisIntegral(energyIndex, csIndex, efficiency);
    }
  }
  return result;
}

/**
 * Evaluating integral with a basis interpolation function
 @param csIndex an index of a corresponding cross section point
*/
double Interpolator::evalIntegralBasis(int csIndex) const {
  double result = 0;
  /**
   * Loop over all interpolation ranges
   */
  for (const auto& rinterp : _rangeInterpolators) {
    if(rinterp.get()->hasCSIndex(csIndex)) {
      /**
       * Evaluate contribution of each interpolation range
       */
      result += rinterp.get()->evalIntegralBasis(csIndex);
    }
  }
  return result;
}
