#include <iostream>
#include <fstream>
#include "Interpolator.hpp"

Interpolator::Interpolator() {}
Interpolator::Interpolator(const Interpolator& interp):
    _rangeInterpolators(interp._rangeInterpolators) {}

Interpolator::Interpolator(const Eigen::VectorXd& cmEnergies,
                           double thresholdEnergy):
    Interpolator(defaultInterpRangeSettings(cmEnergies.rows()),
                 cmEnergies, thresholdEnergy) {}

// TO DO: check ranges and add corresponding exceptions
Interpolator::Interpolator(
    const std::vector<std::tuple<bool, int, int>>& interpRangeSettings,
    const Eigen::VectorXd& cmEnergies, double thresholdEnergy) {
  Eigen::VectorXd extCMEnergies(cmEnergies.rows() + 1);
  extCMEnergies.tail(cmEnergies.rows()) = cmEnergies;
  extCMEnergies(0) = thresholdEnergy;
  _rangeInterpolators.reserve(interpRangeSettings.size());
  for (const auto& el : interpRangeSettings) {
    _rangeInterpolators.push_back(RangeInterpolator(el, extCMEnergies));
  }
}

Interpolator::Interpolator(const std::string& pathToJSON,
                           const Eigen::VectorXd& cmEnergies,
                           double thresholdEnergy):
    Interpolator(Interpolator::loadInterpRangeSettings(pathToJSON),
                 cmEnergies, thresholdEnergy) {}

std::vector<std::tuple<bool, int, int>> Interpolator::loadInterpRangeSettings(
    const std::string& pathToJSON) {
  std::ifstream fl(pathToJSON);
  json s;
  fl >> s;
  std::vector<std::tuple<bool, int, int>> result;
  result.reserve(s.size());
  for (const auto& el : s) {
    result.push_back(std::make_tuple( el[0].get<bool>(),
                                      el[1].get<int>(),
                                      el[2].get<int>()));
  }
  fl.close();
  return result;
}

Interpolator::~Interpolator() {}

double Interpolator::basisEval(int csIndex, double energy) const {
  if (energy <= getMinEnergy()) {
    return 0;
  }
  if (energy > getMaxEnergy()) {
    return _rangeInterpolators.back().basisEval(csIndex, getMaxEnergy());
  }
  double result = 0;
  for (const auto& rinterp : _rangeInterpolators) {
    if (rinterp.hasCSIndex(csIndex) && rinterp.isEnergyInRange(energy)) {
      result += rinterp.basisEval(csIndex, energy);
    }
  }
  return result;
}

double Interpolator::basisDerivEval(int csIndex, double energy) const {
  if (energy <= getMinEnergy()) {
    return 0;
  }
  if (energy > getMaxEnergy()) {
    return _rangeInterpolators.back().basisDerivEval(csIndex, getMaxEnergy());
  }
  double result = 0;
  for (const auto& rinterp : _rangeInterpolators) {
    if (rinterp.hasCSIndex(csIndex) && rinterp.isEnergyInRange(energy)) {
      result += rinterp.basisDerivEval(csIndex, getMaxEnergy());
    }
  }
  return result;
}

double Interpolator::eval(const Eigen::VectorXd& y, double energy) const {
  if (energy <= getMinEnergy()) {
    return 0;
  }
  if (energy > getMaxEnergy()) {
    return _rangeInterpolators.back().eval(y, getMaxEnergy());
  }
  double result = 0;
  for (const auto& rinterp : _rangeInterpolators) {
    if (rinterp.isEnergyInRange(energy)) {
      result += rinterp.eval(y, energy);
    }
  }
  return result;
}

double Interpolator::derivEval(const Eigen::VectorXd& y, double energy) const {
  if (energy <= getMinEnergy()) {
    return 0;
  }
  if (energy > getMaxEnergy()) {
    return _rangeInterpolators.back().derivEval(y, getMaxEnergy());
  }
  double result = 0;
  for (const auto& rinterp : _rangeInterpolators) {
    if (rinterp.isEnergyInRange(energy)) {
      result += rinterp.derivEval(y, energy);
    }
  }
  return result;
}

std::vector<std::tuple<bool, int, int>>
Interpolator::defaultInterpRangeSettings(std::size_t n) {
  return std::vector<std::tuple<bool, int, int>>(
      1, std::tuple<bool, int, int>(false, 0, n));
}

double Interpolator::getMinEnergy() const {
  return _rangeInterpolators.begin()->getMinEnergy();
}

double Interpolator::getMaxEnergy() const {
  return _rangeInterpolators.back().getMaxEnergy();
}
