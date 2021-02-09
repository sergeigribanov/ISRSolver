#ifndef __RANGE_ITERPOLATOR_HPP__
#define __RANGE_ITERPOLATOR_HPP__
#include <gsl/gsl_spline.h>
#include <functional>
#include <memory>
#include <tuple>
#include <vector>
#include <Eigen/Dense>

class RangeInterpolator {
 public:
  RangeInterpolator(const RangeInterpolator&);
  RangeInterpolator(const std::tuple<bool, int, int>&,
                    const Eigen::VectorXd&);
  virtual ~RangeInterpolator();
  double evalKuraevFadinBasisIntegral(int, int, const std::function<double(double, double)>&) const;
  double basisEval(int, double) const;
  double basisDerivEval(int, double) const;
  double eval(const Eigen::VectorXd&, double) const;
  double derivEval(const Eigen::VectorXd&, double) const;
  int getBeginIndex() const;
  double getMinEnergy() const;
  double getMaxEnergy() const;
  bool hasCSIndex(int) const;
  bool isEnergyInRange(double) const;
 private:
  double _evalKuraevFadinIntegralCSpline(
      int, int, const std::function<double(double, double)>&) const;
  double _evalKuraevFadinIntegralLinear(
      int, int, const std::function<double(double, double)>&) const;
  bool _cspline;
  int _beginIndex;
  double _minEnergy;
  double _maxEnergy;
  std::vector<std::shared_ptr<gsl_interp_accel>> _acc;
  std::vector<std::shared_ptr<gsl_spline>> _spline;

};

#endif
