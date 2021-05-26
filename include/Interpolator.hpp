#ifndef _INTERPOLATOR_HPP_
#define _INTERPOLATOR_HPP_
#include <string>
#include <exception>
#include <memory>
#include <nlohmann/json.hpp>
#include "BaseRangeInterpolator.hpp"

using json = nlohmann::json;

/**
 * Exception: wrong interpolation range
 */
typedef struct : std::exception {
  const char* what() const noexcept {
    return "[!] Wrong interpolation range.\n";
  }
} InterpRangeException;

/**
 * This class is designed in order to interpolate unknown Born
 * cross section points
 */
class Interpolator {
 public:
  /**
   * Default constructor
   */
  Interpolator();
  /**
   * Copy constructor
   */
  Interpolator(const Interpolator&);
  /**
   * Constructor
   * @param cmEnergies center-of-mass energies
   * @param thresholdEnergy a threshold energy
   */
  Interpolator(
      const Eigen::VectorXd& cmEnergies,
      double thresholdEnergy);
  /**
   * Constructor
   * @param interpRangeSettings an interpolation settings
   * @param cmEnergies center-of-mass energies
   * @param thresholdEnergy a threshold energy
   */
  Interpolator(
      const std::vector<std::tuple<bool, int, int>>&
      interpRangeSettings,
      const Eigen::VectorXd& cmEnergies,
      double thresholdEnergy) noexcept(false);
  /**
   * Constructor
   * @param pathToJSON a path to .json file with interpolation settings
   * @param cmEnergies a vector of center-of-mass energies
   * @param thresholdEnergy a threshold energy
   */
  Interpolator
  (const std::string&  pathToJSON,
   const Eigen::VectorXd& cmEnergies,
   double thresholdEnergy);

  Interpolator(const json& obj,
               const Eigen::VectorXd& cmEnergies,
               double thresholdEnergy);

  /**
   * Destructor
   */
  virtual ~Interpolator();
  /**
   * Kuraev-Fadin convolution with a basis interpolation function
   * that corresponds to the cross section point indx csIndex.
   * @param energyIndex an index of center-of-mass energy at
   * which convolution is evaluated
   * @param csIndex an index of cross-section point that corresponds
   * to considered basis interpolation function
   * @param efficiency a detection efficiency
   */
  double evalKuraevFadinBasisIntegral(
      int energyIndex,
      int csIndex,
      const std::function<double(double, double)>& efficiency) const;
  double evalBasisSConvolution(
      int csIndex,
      const std::function<double(double)>& convKernel) const;
  /**
   * Evaluating integral with a basis interpolation function
   @param csIndex an index of a corresponding cross section point
   */
  double evalIntegralBasis(int csIndex) const;
  /**
   * Evaluating basis interpolation function
   */
  double basisEval(int csIndex, double energy) const;
  /**
   * Evaluating derivative of a basis interpolation function
   * @param csIndex an index of a corresponding cross section
   point
   @param energy is a center-of-mass energy
   */
  double basisDerivEval(int csIndex, double energy) const;
  /**
   * Evaluating interpolation y data value at a certain
   center-of-mass energy
   * @param y a data
   * @param energy a center-of-mass energy
   */
  double eval(const Eigen::VectorXd& y, double energy) const;
  /**
   * Evaluating interpolation derivative
   * @param y a data
   * @param energy a center-of-mass energy
   */
  double derivEval(
      const Eigen::VectorXd& y,
      double energy) const;
  /**
   * Getting the minimum center-of-mass energy associated with
   * the considered interpolator
   */
  double getMinEnergy() const;
  /**
   * Getting the maximum center-of-mass energy associated with
   * the considered interpolator
   */
  double getMaxEnergy() const;
  /**
   * Getting the minimum center-of-mass energy
   * that belongs interpolation range, which contains
   * the center-of-mass energy segment with index csIndex
   * @param csIndex an index of center-of-mass energy segment
  */
  double getMinEnergy(int csIndex) const;
  /**
   * Getting the maximum center-of-mass energy
   * that belongs interpolation range, which contains
   * the center-of-mass energy segment with index csIndex
   * @param csIndex an index of center-of-mass energy segment
   */
  double getMaxEnergy(int csIndex) const;

  static std::vector<std::tuple<bool, int, int>> loadInterpRangeSettingJSON(const json& obj);
  /**
   * Load interpolation settings from file
   * @param pathToJSON a path to a .json file with settings
   */
  static std::vector<std::tuple<bool, int, int>>
  loadInterpRangeSettings(
      const std::string& pathToJSON);
  /**
   * Use default interpolation settings (piecewise linear interpolation)
   * @param n a number of cross section points
   */
  static std::vector<std::tuple<bool, int, int>>
  defaultInterpRangeSettings(int n);
 private:
  /**
   * Check validity of range interpolation settings
   * @param sortedInterpRangeSettings an interpolation settings
   * @param cmEnergies a vector of center-of-mass energies
   */
  static bool checkRangeInterpSettings(
      const std::vector<std::tuple<bool, int, int>>&
      sortedInterpRangeSettings,
      const Eigen::VectorXd& cmEnergies);
  /**
   * Interpolation settings
   */
  std::vector<std::shared_ptr<BaseRangeInterpolator>> _rangeInterpolators;
};

#endif
