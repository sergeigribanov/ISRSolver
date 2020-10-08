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
  void solveTikhonov();
  void save(const std::string& outputPath, const OutputOptions& outputOpts);
  void setInterpSettings(const std::vector<InterpPtSettings>&) noexcept(false);
  void setInterpSettings(const std::string&) noexcept(false);

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
  double evalRegFuncNorm2(const Eigen::VectorXd&) const;
  double evalDifferenceNorm2() const;
  double evalSolDSolNorm2() const;
  double evalSolNorm2() const;
  double evalDSolNorm2() const;
  Eigen::VectorXd evalRegFuncGradNorm2(const Eigen::VectorXd&) const;
  Eigen::MatrixXd evalHessian() const;
  TF1* createInterpFunction() const;
  TF1* createDerivativeInterpFunction(std::size_t, const std::string&) const;
  void setDefaultInterpSettings();
  Eigen::RowVectorXd interpProjector(double) const;
  Eigen::RowVectorXd interpDerivativeProjector(double, std::size_t) const;
  Eigen::MatrixXd interpPointWiseDerivativeProjector() const;
  Eigen::RowVectorXd polIntegralOp(int) const;
  Eigen::RowVectorXd evalPolA(int) const;
  Eigen::MatrixXd evalIntegralMatrix() const;
  Eigen::RowVectorXd evalScalarProductOperator() const;
  friend double objectiveFCN(unsigned, const double*, double*, void*);
  double _alpha;
  double _beta;
  
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
