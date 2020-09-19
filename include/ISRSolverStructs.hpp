#ifndef __ISRSOLVER_STRUCTS_HPP__
#define __ISRSOLVER_STRUCTS_HPP__

#include <exception>

typedef struct {
  double cmEnergy;
  double cs;
  double cmEnergyError;
  double csError;
} CSData;

typedef struct {
  Eigen::VectorXd cmEnergy;
  Eigen::VectorXd s;
  Eigen::VectorXd cs;
  Eigen::VectorXd cmEnergyError;
  Eigen::MatrixXd invCSErrMatrix;
} CSVecData;

typedef struct {
  std::size_t numCoeffs;
  std::size_t nPointsLeft;
} InterpPtSettings;

typedef struct {
  std::string measuredCSGraphName;
  double thresholdEnergy;
  bool energyUnitMeVs;
} InputOptions;

typedef struct {
  std::string measuredCSGraphName;
  std::string bornCSGraphName;
} OutputOptions;

typedef struct : std::exception {
  const char* what() const noexcept {
    return "[!] Wrong size of interpolation settings vector.\n";
  }
} InterpSettingsSizeException;

#endif
