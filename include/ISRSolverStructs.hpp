#ifndef _ISRSOLVER_STRUCTS_HPP_
#define _ISRSOLVER_STRUCTS_HPP_

#include <exception>

/**
 * Cross section data at a single point
 */
typedef struct {
  /**
   * Center-of-mass energy
   */
  double cmEnergy;
  /**
   * Cross section
   */
  double cs;
  /**
   * Center-of-mass energy error
   */
  double cmEnergyError;
  /**
   * Cross section error
   */
  double csError;
} CSData;

/**
 * Cross section at all points
 */
typedef struct {
  /**
   * Center-of-mass energy values
   */
  Eigen::VectorXd cmEnergy;
  /**
   * Center-of-mass energy error values
   */
  Eigen::VectorXd cmEnergyError;
  /**
   * Cross section values
   */
  Eigen::VectorXd cs;
  /**
   * Cross section error values
   */
  Eigen::VectorXd csError;
} CSVecData;

/**
 * Input options
 */
typedef struct {
  /**
   * A name of a detection efficiency object (TEfficiency)
   */
  std::string efficiencyName;
  /**
   * A name of a visible cross section object (TGraphErrors)
   */
  std::string visibleCSGraphName;
  /**
   * A threshold energy
   */
  double thresholdEnergy;
} InputOptions;

/**
 * Output options
 */
typedef struct {
  /**
   * A name of a visible cross section object (TGraphErrors)
   */
  std::string visibleCSGraphName;
  /**
   * A name a Born cross section object (TGraphErrors)
   */
  std::string bornCSGraphName;
} OutputOptions;

/**
 * Interpolation settings exception
 */
typedef struct : std::exception {
  const char* what() const noexcept {
    return "[!] Wrong size of interpolation settings vector.\n";
  }
} InterpSettingsSizeException;

#endif
