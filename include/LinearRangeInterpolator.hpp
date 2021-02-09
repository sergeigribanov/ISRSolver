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
  private:
  double _c00(int) const;
  double _c01(int) const;
  double _c10(int) const;
  double _c11(int) const;
  Eigen::VectorXd _extCMEnergies;
};

#endif
