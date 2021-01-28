#ifndef __INTERPOLATOR_HPP__
#define __INTERPOLATOR_HPP__
#include <string>
#include <nlohmann/json.hpp>
#include "RangeInterpolator.hpp"

using json = nlohmann::json;

class Interpolator {
 public:
  Interpolator();
  Interpolator(const Interpolator&);
  Interpolator(const Eigen::VectorXd&, double);
  Interpolator(const std::vector<std::tuple<bool, int, int>>&,
               const Eigen::VectorXd&, double);
  Interpolator(const std::string&, const Eigen::VectorXd&, double);
  virtual ~Interpolator();
  double basisEval(int, double) const;
  double basisDerivEval(int, double) const;
  double eval(const Eigen::VectorXd&, double) const;
  double derivEval(const Eigen::VectorXd&, double) const;
  double getMinEnergy() const;
  double getMaxEnergy() const;
  static std::vector<std::tuple<bool, int, int>>
  loadInterpRangeSettings(const std::string&);
  static std::vector<std::tuple<bool, int, int>>
  defaultInterpRangeSettings(int);
 private:
  std::vector<RangeInterpolator> _rangeInterpolators;
};

#endif
