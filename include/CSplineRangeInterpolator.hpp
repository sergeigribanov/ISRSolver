#ifndef __CSPLINE_RANGE_INTERPOLATOR_HPP__
#define __CSPLINE_RANGE_INTERPOLATOR_HPP__
#include <vector>
#include <memory>
#include <gsl/gsl_spline.h>
#include "BaseRangeInterpolator.hpp"

class CSplineRangeInterpolator : public BaseRangeInterpolator {
 public:
  CSplineRangeInterpolator(int, int, const Eigen::VectorXd&);
  CSplineRangeInterpolator(const CSplineRangeInterpolator&);
  virtual ~CSplineRangeInterpolator();
  double basisEval(int, double) const override final;
  double basisDerivEval(int, double) const override final;
  double evalKuraevFadinBasisIntegral(
      int, int, const std::function<double(double, double)>&) const override final;
  double evalIntegralBasis(int) const override final;
 private:
  std::vector<std::shared_ptr<gsl_interp_accel>> _acc;
  std::vector<std::shared_ptr<gsl_spline>> _spline;
};

#endif
