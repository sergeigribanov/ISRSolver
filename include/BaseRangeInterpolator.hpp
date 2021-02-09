#ifndef __BASE_RANGE_ITERPOLATOR_HPP__
#define __BASE_RANGE_ITERPOLATOR_HPP__
#include <functional>
#include <Eigen/Dense>

class BaseRangeInterpolator {
 public:
  BaseRangeInterpolator(int, int, const Eigen::VectorXd&);
  BaseRangeInterpolator(const BaseRangeInterpolator&);
  bool hasCSIndex(int) const;
  bool isEnergyInRange(double) const;
  int getBeginIndex() const;
  int getNumberOfSegments() const;
  double eval(const Eigen::VectorXd&, double) const;
  double derivEval(const Eigen::VectorXd&, double) const;
  double getMinEnergy() const;
  double getMaxEnergy() const;
  virtual double basisEval(int, double) const = 0;
  virtual double basisDerivEval(int, double) const = 0;
  virtual double evalKuraevFadinBasisIntegral(
      int, int, const std::function<double(double, double)>&) const = 0;
 protected:
  static int _evalNumOfSegments(int, int);
  static int _evalBeginIndex(int, int);
  int _numberOfSegments;
  int _beginIndex;
  double _minEnergy;
  double _maxEnergy;
};

#endif
