#ifndef _ISRSOLVER_STRUCTS_HPP_
#define _ISRSOLVER_STRUCTS_HPP_

#include <exception>
#include <Eigen/Dense>

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
 * Chi-square model test arguments
 * @see chi2TestModel
 */
typedef struct {
  /**
   * A number of numerical experiments
   */
  int n;
  /**
   * An initial value of chi-square fit parameter
   */
  double initialChi2Ampl;
  /**
   * A path to model .root file
   */
  std::string modelPath;
  /**
   * A name of model visible cross section
   object (TGraphErrors)
   */
  std::string modelVCSName;
  /**
   * A name of model Born cross section
   object (TGraphErrors)
   */
  std::string modelBCSName;
  /**
   * A path to .root file with results
   */
  std::string outputPath;
} Chi2TestModelArgs;


/**
 * Chi-square test arguments
 * @see chi2TestData
 * @see chi2Test
 */
typedef struct {
  /**
   * A number of numerical experiments
   */
  int n;
  /**
   * An initial value of chi-square fit parameter
   */
  double initialChi2Ampl;
  /**
   * A path to .root file with results
   */
  std::string outputPath;
} Chi2TestArgs;

/**
 * Ratio model test arguments
 * (Ratio of a numerical solution to a model Born cross section)
 * @see ratioTestModel
 */
typedef struct {
  /**
   * A number of numerical experiments
   */
  int n;
  /**
   * A path to model .root file
   */
  std::string modelPath;
  /**
   * A name of model visible cross section
   object (TGraphErrors)
  */
  std::string modelVCSName;
  /**
   * A name of model Born cross section
   object (TGraphErrors)
  */
  std::string modelBCSName;
  /**
   * A path to .root file with results
   */
  std::string outputPath;
} RatioTestModelArgs;

/**
 * Interpolation settings exception
 */
typedef struct : std::exception {
  const char* what() const noexcept {
    return "[!] Wrong size of interpolation settings vector.\n";
  }
} InterpSettingsSizeException;

#endif
