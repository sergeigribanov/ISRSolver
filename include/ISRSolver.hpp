#ifndef __ISRSOLVER_HPP__
#define __ISRSOLVER_HPP__
#include <TF1.h>
#include <TGraphErrors.h>
#include <TMatrixT.h>

#include <Eigen/Dense>
#include <string>
#include <utility>
#include <vector>

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

class ISRSolver {
 public:
  ISRSolver(const std::string& inputPath,
	    const InputOptions& inputOpts);
  virtual ~ISRSolver();
  std::size_t getN() const;
  // void solve();
  // void save(const std::string& outputPath,
  // 	    const OutputOptions& outputOpts);
  void testPrint() const;
 private:
  double getXmin(int) const;
  double getXmax(int) const;
  double getNumCoeffs(int) const;
  double getNPointsLeft(int) const;
  Eigen::MatrixXd permutation(int) const;
  Eigen::MatrixXd interpInvMatrix(int) const;
  Eigen::MatrixXd polConvKuraevFadinMatrix(int) const;
  Eigen::MatrixXd evalA(int) const;
  Eigen::MatrixXd evalEqMatrix() const;
  void setDefaultInterpSettings();
  InputOptions _inputOpts;
  double _sT;
  std::size_t _n;
  std::vector<InterpPtSettings> _interpSettings;
  CSVecData _measuredCSData;
  std::vector<double> _xPoints;
};
#endif
