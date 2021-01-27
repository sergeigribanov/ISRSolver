#include <fstream>
#include "Interpolator.hpp"

Interpolator::Interpolator() {}
Interpolator::Interpolator(const Interpolator& interp):
    _rangeInterpolators(interp._rangeInterpolators) {}


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
    Interpolator(Interpolator::loadInterpSettings(pathToJSON),
                 cmEnergies, thresholdEnergy) {}

std::vector<std::tuple<bool, int, int>> Interpolator::loadInterpSettings(
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
  double result = 0;
  for (const auto& rinterp : _rangeInterpolators) {
    if (rinterp.hasCSIndex(csIndex) && rinterp.isEnergyInRange(energy)) {
      result += rinterp.basisEval(csIndex, energy);
    }
  }
  return result;
}

double Interpolator::basisDerivEval(int csIndex, double energy) const {
  double result = 0;
  for (const auto& rinterp : _rangeInterpolators) {
    if (rinterp.hasCSIndex(csIndex) && rinterp.isEnergyInRange(energy)) {
      result += rinterp.basisDerivEval(csIndex, energy);
    }
  }
  return result;
}

double Interpolator::eval(const Eigen::VectorXd& y, double energy) const {
  double result = 0;
  for (const auto& rinterp : _rangeInterpolators) {
    if (rinterp.isEnergyInRange(energy)) {
      result += rinterp.eval(y, energy);
    }
  }
  return result;
}

double Interpolator::derivEval(const Eigen::VectorXd& y, double energy) const {
  double result = 0;
  for (const auto& rinterp : _rangeInterpolators) {
    if (rinterp.isEnergyInRange(energy)) {
      result += rinterp.derivEval(y, energy);
    }
  }
  return result;
}
