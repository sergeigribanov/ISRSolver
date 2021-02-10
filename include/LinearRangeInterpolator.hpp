#ifndef __LINEAR_RANGE_INTERPOLATOR_HPP__
#define __LINEAR_RANGE_INTERPOLATOR_HPP__
#include "BaseRangeInterpolator.hpp"

class LinearRangeInterpolator : public BaseRangeInterpolator {
 public:
  LinearRangeInterpolator(int, int, const Eigen::VectorXd&);
  LinearRangeInterpolator(const LinearRangeInterpolator&);
  virtual ~LinearRangeInterpolator();
  double basisEval(int, double) const override final;
  double basisDerivEval(int, double) const override final;
  double evalKuraevFadinBasisIntegral(
      int, int, const std::function<double(double, double)>&) const override final;
  double _evalKuraevFadinBasisIntegralFirstTriangle(
      int, int, const std::function<double(double, double)>&) const;
  double _evalKuraevFadinBasisIntegralSecondTriangle(
      int, int, const std::function<double(double, double)>&) const;
  double evalIntegralBasis(int) const override final;
  private:
  double _evalIntegralBasisFirstTriangle(int) const;
  double _evalIntegralBasisSecondTriangle(int) const;
  double _c00(int) const;
  double _c01(int) const;
  double _c10(int) const;
  double _c11(int) const;
  Eigen::VectorXd _extCMEnergies;
  static std::function<double(double)> _fcn0;
  static std::function<double(double)> _fcn1;
};

#endif
