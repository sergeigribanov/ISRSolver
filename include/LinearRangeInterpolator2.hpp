#ifndef _LINEAR_RANGE_INTERPOLATOR2_HPP_
#define _LINEAR_RANGE_INTERPOLATOR2_HPP_
#include "BaseRangeInterpolator2.hpp"

/**
 * This class is designed in order to perform a linear spline
 interpolation of the numerical solution in a certain range
 */
class LinearRangeInterpolator2 : public BaseRangeInterpolator2 {
 public:
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
  LinearRangeInterpolator2(
      int rangeIndexMin,
      int rangeIndexMax,
      const Eigen::VectorXd& extCMEnergies);
  /**
   * Copy constructor
   */
  LinearRangeInterpolator2(const LinearRangeInterpolator2&);
  /**
   * Destructor
   */
  virtual ~LinearRangeInterpolator2();
  /**
   * Evaluate basis interpolation around center-of-mass energy that
     corresponds to the cross section point with index the csIndex
     * @param csIndex a cross section point index
     * @param energy a center-of-mass energy at which interpolation
     is evaluated
  */
  double basisEval(int csIndex, double energy) const override final;
  /**
   * Evaluate basis interpolation derivative around center-of-mass energy that
   corresponds to the cross section point with the csIndex
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
      const std::function<double(double, std::size_t)>& efficiency) const override final;
  virtual double evalBasisSConvolution(
      int csIndex,
      const std::function<double(double)>& convKernel) const override final;
  virtual double evalBasisSConvolution(
      int csIndex,
      const std::function<double(double)>& convKernel,
      double s_min, double s_max) const override final;
  /**
   * Evaluate Kuraev-Fadin convolution inside the first triangle
   with basis interpolation
   * @param energyIndex a center-of-mass energy index
   * @param csIndex a cross section point index
   * @param efficiency a detection efficiency
   */
  double _evalKuraevFadinBasisIntegralFirstTriangle(
      int energyIndex, int csIndex,
      const std::function<double(double, std::size_t)>& efficiency) const;
  /**
   * Evaluate Kuraev-Fadin convolution inside the second triangle
   with basis interpolation
   * @param energyIndex a center-of-mass energy index
   * @param csIndex a cross section point index
   * @param efficiency a detection efficiency
   */
  double _evalKuraevFadinBasisIntegralSecondTriangle(
      int energyIndex, int csIndex,
      const std::function<double(double, std::size_t)>& efficiency) const;
  double _evalBasisSConvolutionFirstTriangle(
      int csIndex,
      const std::function<double(double)>& convKernel) const;
  double _evalBasisSConvolutionSecondTriangle(
      int csIndex,
      const std::function<double(double)>& convKernel) const;
  double _evalBasisSConvolutionFirstTriangle(
      int csIndex,
      const std::function<double(double)>& convKernel,
      double s_min, double s_max) const;
  double _evalBasisSConvolutionSecondTriangle(
      int csIndex,
      const std::function<double(double)>& convKernel,
      double s_min, double s_max) const;
  /** Evaluate integral of basis interpolation function that correspond to
      a csIndex-th cross section point
      * @param csIndex a cross section point index
      */
  double evalIntegralBasis(int csIndex) const override final;
  private:
  /**
   * Integral of a basis interpolation function inside the first
   triangle
   *
   */
  double _evalIntegralBasisFirstTriangle(int csIndex) const;
  /**
   * Evaluate integral of a basis interpolation function inside
   the second triangle
   */
  double _evalIntegralBasisSecondTriangle(int csIndex) const;
  /**
   * Auxiliary coefficients:
   _c00, _c01, _c10, _c11
   *
   */
  double _c00(int csIndex) const;
  double _c01(int csIndex) const;
  double _c10(int csIndex) const;
  double _c11(int csIndex) const;
  Eigen::VectorXd _extCMEnergies;
  /**
   * f(x) = 1
   */
  static std::function<double(double)> _fcn0;
  /**
   * f(x) = x
   */
  static std::function<double(double)> _fcn1;
};

#endif
