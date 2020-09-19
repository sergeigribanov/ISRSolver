#ifndef __ISRSOLVER_HPP__
#define __ISRSOLVER_HPP__

#include <TF1.h>

#include <Eigen/Dense>
#include <string>
#include <vector>

#include "ISRSolverStructs.hpp"

class ISRSolver {
 public:
  ISRSolver(const std::string& inputPath, const InputOptions& inputOpts);
  virtual ~ISRSolver();
  std::size_t getN() const;
  void solve();
  void save(const std::string& outputPath, const OutputOptions& outputOpts);
  void setInterpSettings(const std::vector<InterpPtSettings>&) noexcept(false);

 private:
  double getXmin(int, int) const;
  double getXmax(int, int) const;
  double getNumCoeffs(int) const;
  double getNPointsLeft(int) const;
  Eigen::MatrixXd permutation(int) const;
  Eigen::MatrixXd interpInvMatrix(int) const;
  Eigen::MatrixXd polConvKuraevFadinMatrix(int) const;
  Eigen::MatrixXd evalA(int) const;
  Eigen::MatrixXd evalEqMatrix() const;
  TF1* createInterpFunction() const;
  void setDefaultInterpSettings();
  InputOptions _inputOpts;
  double _sT;
  std::size_t _n;
  std::vector<InterpPtSettings> _interpSettings;
  CSVecData _measuredCSData;
  Eigen::MatrixXd _integralOperatorMatrix;
  Eigen::VectorXd _bornCS;
  Eigen::MatrixXd _invBornCSErrorMatrix;
};
#endif
