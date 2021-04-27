#ifndef _CSPLINE_RANGE_INTERPOLATOR_HPP_
#define _CSPLINE_RANGE_INTERPOLATOR_HPP_
#include <vector>
#include <memory>
#include <gsl/gsl_spline.h>
#include "BaseRangeInterpolator.hpp"

/*
 * This class is designed in order to perform a cubic spline interpolation
 * of the numerical solution in a certain range
 */
class CSplineRangeInterpolator : public BaseRangeInterpolator {
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
  CSplineRangeInterpolator(int rangeIndexMin,
                           int rangeIndexMax,
                           const Eigen::VectorXd& extCMEnergies);
  /**
   * Copy constructor
   */
  CSplineRangeInterpolator(const CSplineRangeInterpolator&);
  /**
   * Destructor
   */
  virtual ~CSplineRangeInterpolator();
  /**
   * Evaluate basis interpolation around center-of-mass energy that
   * corresponds to the cross section point with index the csIndex
   * @param csIndex a cross section point index
   * @param energy a center-of-mass energy at which interpolation
   * is evaluated
   */
  double basisEval(int csIndex, double energy) const override final;
  /**
   * Evaluate basis interpolation derivative around center-of-mass energy that
   * corresponds to the cross section point with the csIndex
   * @param csIndex a cross section point index
   * @param energy a center-of-mass energy a which interpolation is evaluated
   */
  double basisDerivEval(int csIndex, double energy) const override final;
  /**
   * Evaluate Kuraev-Fadin convolution with basis interpolation
   * @param energyIndex a center-of-mass energy index
   * @param csIndex a cross section point index
   * @param efficiency a detection efficiency
   */
  double evalKuraevFadinBasisIntegral(
      int energyIndex,
      int csIndex,
      const std::function<double(double, double)>& efficiency) const override final;
  /** Evaluate integral of basis interpolation function that correspond to
   * a csIndex-th cross section point
   * @param csIndex a cross section point index
   */
  double evalIntegralBasis(int csIndex) const override final;
 private:
  /**
   * GSL accelerators
   */
  std::vector<std::shared_ptr<gsl_interp_accel>> _acc;
  /**
   * GSL splines
   */
  std::vector<std::shared_ptr<gsl_spline>> _spline;
};

#endif
