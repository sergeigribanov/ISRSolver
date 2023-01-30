#ifndef _BASE_RANGE_ITERPOLATOR2_HPP_
#define _BASE_RANGE_ITERPOLATOR2_HPP_
#include <functional>
#include <Eigen/Dense>

/**
 * Base class of a range interpolator
 */
class BaseRangeInterpolator2 {
 public:
  /**
   * Constructor
   * @param rangeIndexMin a minimum index of sub range (a range between
   * two neighboring center-of-mass energy points) that used in interpolation
   * range
   * @param rangeIndexMax a maximum index of sub range that used in
   * interpolation range
   * @param extCMEnergies extended vector of center-of-mass energies (contains
   * threshold energy)
   */
  BaseRangeInterpolator2(
      int rangeIndexMin,
      int rangeIndexMax,
      const Eigen::VectorXd& extCMEnergies);
  /**
   * Copy constructor
   */
  BaseRangeInterpolator2(const BaseRangeInterpolator2&);
  /**
   * This method returns true if cross section point
   * with csIndex belongs this interpolation sub range
   * @param csIndex an index of a cross section point
   */
  bool hasCSIndex(int csIndex) const;
  /**
   * This method returns true if energy belongs
   * to the considered interpolation sub range
   * @param energy a center-of-mass energy
   */
  bool isEnergyInRange(double energy) const;
  /**
   * This method returns the index of the first segment
   * (a segment with lowest center-of-mass energies) that belongs
   * to the considered sub range
   */
  int getBeginIndex() const;
  /**
   * This method returns a number of center-of-mass energy
   * segments (from one
   * cross section point to another) in the considered
   * sub range.
   * The segment indices are defined as follows.
   * If index=0, the center-of-mass energy segment
   * between threshold energy and the first center-of-mass
   * energy is considered. In the case when index=1,
   * the center-of-mass energy segment between the first
   * and the second center-of-mass energy points is considered.
   */
  int getNumberOfSegments() const;
  /**
   * This method is used to eval interpolation
   * value at the certain center-of-mass energy
   * @param y a data
   * @param energy a center-of-mass energy
   */
  double eval(const Eigen::VectorXd& y,
              double energy) const;
  /**
   * This method is used to eval interpolation
   * derivative value at the certain center-of-mass energy
   * @param y a data
   * @param energy a center-of-mass energy
   */
  double derivEval(const Eigen::VectorXd& y,
                   double energy) const;
  /**
   * This method is used to get minimum energy
   */
  double getMinEnergy() const;
  /**
   * This method is used to get maximum energy
   */
  double getMaxEnergy() const;
  /**
   * Evaluate basis interpolation around center-of-mass energy that
   * corresponds to the cross section point with index the csIndex
   * @param csIndex a cross section point index
   * @param energy a center-of-mass energy at which interpolation
   * is evaluated
  */
  virtual double basisEval(int csIndex, double energy) const = 0;
  /**
   * Evaluate basis interpolation derivative around center-of-mass energy that
   * corresponds to the cross section point with the csIndex
   * @param csIndex a cross section point index
   * @param energy a center-of-mass energy a which interpolation is evaluated
   */
  virtual double basisDerivEval(int csIndex, double energy) const = 0;
  /**
   * Evaluate Kuraev-Fadin convolution with basis interpolation
   * @param energyIndex a center-of-mass energy index
   * @param csIndex a cross section point index
   * @param efficiency a detection efficiency
   */
  virtual double evalKuraevFadinBasisIntegral(
      int energyIndex,
      int csIndex,
      const std::function<double(double, std::size_t)>& efficiency) const = 0;
  virtual double evalBasisSConvolution(
      int csIndex,
      const std::function<double(double)>& convKernel) const = 0;
  virtual double evalBasisSConvolution(
      int csIndex,
      const std::function<double(double)>& convKernel,
      double s_min, double s_max) const = 0;
  /** Evaluate integral of basis interpolation function that correspond to
   * a csIndex-th cross section point
   * @param csIndex a cross section point index
   */
  virtual double evalIntegralBasis(int csIndex) const = 0;
 protected:
  /**
   * Evaluating number of segments of the current
   * interpolation sub range
   * @param rangeIndexMin a minimum sub range index
   * @param rangeIndexMax a maximum sub range index
   */
  static int _evalNumOfSegments(
      int rangeIndexMin,
      int rangeIndexMax);
  /**
   * Evaluating the begin index
   * @param rangeIndexMin a minimum sub range index
   * @param rangeIndexMax a maximum sub range index
   */
  static int _evalBeginIndex(
      int rangeIndexMin,
      int rangeIndexMax);
  /**
   * Number of sub ranges
   */
  int _numberOfSegments;
  /**
   * The lowest segment index in the
   * considered interpolation sub range.
   */
  int _beginIndex;
  /**
   * Minimum energy
   */
  double _minEnergy;
  /**
   * Maximum energy
   */
  double _maxEnergy;
};

#endif
