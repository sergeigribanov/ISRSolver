#ifndef __RANGE_ITERPOLATOR_HPP__
#define __RANGE_ITERPOLATOR_HPP__
#include <gsl/gsl_spline.h>
#include <tuple>
#include <vector>
#include <Eigen/Dense>

class RangeInterpolator {
 public:
  RangeInterpolator(const std::tuple<bool, int, int>&,
                    const Eigen::VectorXd&);
  virtual ~RangeInterpolator();
  double basisEval(int, double) const;
  double basisDerivEval(int, double) const;
  double eval(const Eigen::VectorXd&, double) const;
  double derivEval(const Eigen::VectorXd&, double) const;
 private:
  int _beginIndex;
  double _minEnergy;
  double _maxEnergy;
  std::vector<gsl_interp_accel*> _acc;
  std::vector<gsl_spline*> _spline;

};

#endif
