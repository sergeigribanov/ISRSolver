#include <algorithm>
#include <fstream>
#include "Interpolator.hpp"

Interpolator::Interpolator() {}
Interpolator::Interpolator(const Interpolator& interp):
    _rangeInterpolators(interp._rangeInterpolators) {}

Interpolator::Interpolator(const Eigen::VectorXd& cmEnergies,
                           double thresholdEnergy):
    Interpolator(defaultInterpRangeSettings(cmEnergies.rows()),
                 cmEnergies, thresholdEnergy) {}

Interpolator::Interpolator(
    const std::vector<std::tuple<bool, int, int>>& interpRangeSettings,
    const Eigen::VectorXd& cmEnergies, double thresholdEnergy) {
  Eigen::VectorXd extCMEnergies(cmEnergies.rows() + 1);
  extCMEnergies.tail(cmEnergies.rows()) = cmEnergies;
  extCMEnergies(0) = thresholdEnergy;
  auto interpRSsorted = interpRangeSettings;
  std::sort(interpRSsorted.begin(), interpRSsorted.end(),
            [](const std::tuple<bool, int, int>& t1,
               const std::tuple<bool, int, int>& t2) {
              return std::get<1>(t1) < std::get<1>(t2);});
  _rangeInterpolators.reserve(interpRSsorted.size());
  if (checkRangeInterpSettings(interpRSsorted,
                               cmEnergies)) {
    InterpRangeException ex;
    throw ex;
  }
  for (const auto& el : interpRSsorted) {
    _rangeInterpolators.push_back(RangeInterpolator(el, extCMEnergies));
  }
}

bool Interpolator::checkRangeInterpSettings(
    const std::vector<std::tuple<bool, int, int>>&
    sortedInterpRangeSettings,
    const Eigen::VectorXd& cmEnergies) {
  if(std::get<1>(*sortedInterpRangeSettings.begin()) != 0 ||
     std::get<2>(sortedInterpRangeSettings.back()) + 1
     != cmEnergies.rows()) {
    return true;
  }
  int pimax = -1;
  for (const auto& el : sortedInterpRangeSettings) {
    int imin = std::get<1>(el);
    int imax = std::get<2>(el);
    if (imin > imax) {
      return true;
    }
    if (pimax != -1 && pimax + 1 != imin) {
      return true;
    }
    pimax = imax;
  }
  return false;
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
    if (_rangeInterpolators.back().hasCSIndex(csIndex)) {
      return _rangeInterpolators.back().basisEval(csIndex, getMaxEnergy());
    } else {
      return 0;
    }
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
    if (_rangeInterpolators.back().hasCSIndex(csIndex)) {
      return _rangeInterpolators.back().basisDerivEval(csIndex, getMaxEnergy());
    } else {
      return 0;
    }
  }
  double result = 0;
  for (const auto& rinterp : _rangeInterpolators) {
    if (rinterp.hasCSIndex(csIndex) && rinterp.isEnergyInRange(energy)) {
      result += rinterp.basisDerivEval(csIndex, energy);
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
Interpolator::defaultInterpRangeSettings(int n) {
  return std::vector<std::tuple<bool, int, int>>(
      1, std::tuple<bool, int, int>(false, 0, n - 1));
}

double Interpolator::getMinEnergy() const {
  return _rangeInterpolators.begin()->getMinEnergy();
}

double Interpolator::getMaxEnergy() const {
  return _rangeInterpolators.back().getMaxEnergy();
}

double Interpolator::getMinEnergy(int csIndex) const {
  double result = _rangeInterpolators.back().getMinEnergy();
  for (const auto& rinterp : _rangeInterpolators) {
    if(rinterp.hasCSIndex(csIndex)) {
      result = std::min(result, rinterp.getMinEnergy());
    }
  }
  return result;
}

double Interpolator::getMaxEnergy(int csIndex) const {
  double result = _rangeInterpolators[0].getMaxEnergy();
  for (const auto& rinterp : _rangeInterpolators) {
    if(rinterp.hasCSIndex(csIndex)) {
      result = std::max(result, rinterp.getMaxEnergy());
    }
  }
  return result;
}
