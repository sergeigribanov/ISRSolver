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
  Eigen::VectorXd cmEnergyError;
  Eigen::VectorXd cs;
  Eigen::VectorXd csError;
} CSVecData;

typedef struct {
  std::size_t numCoeffs;
  std::size_t nPointsLeft;
} InterpPtSettings;

typedef struct {
  std::string visibleCSGraphName;
  double thresholdEnergy;
  bool energyUnitMeVs;
} InputOptions;

typedef struct {
  std::string visibleCSGraphName;
  std::string bornCSGraphName;
} OutputOptions;

typedef struct : std::exception {
  const char* what() const noexcept {
    return "[!] Wrong size of interpolation settings vector.\n";
  }
} InterpSettingsSizeException;

#endif
