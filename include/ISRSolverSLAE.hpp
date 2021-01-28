#ifndef __ISRSOLVER_SLAE_HPP__
#define __ISRSOLVER_SLAE_HPP__

#include <TF1.h>
#include "Interpolator.hpp"
#include "BaseISRSolver.hpp"

class ISRSolverSLAE : public BaseISRSolver {
 public:
  ISRSolverSLAE(const std::string& inputPath, const InputOptions& inputOpts);
  ISRSolverSLAE(const ISRSolverSLAE&);
  virtual ~ISRSolverSLAE();
  const Eigen::MatrixXd& getIntegralOperatorMatrix() const;
  const Eigen::MatrixXd& getInverseBornCSErrorMatrix() const;
  virtual void solve() override;
  virtual void save(const std::string& outputPath,
                    const OutputOptions& outputOpts) override;
  void setRangeInterpSettings(const std::vector<std::tuple<bool, int, int>>&);
  void setRangeInterpSettings(const std::string&);
  // void setInterpSettings(const std::vector<InterpPtSettings>&) noexcept(false);
  // void setInterpSettings(const std::string&) noexcept(false);

 protected:
  // double _getXmin(int, int) const;
  // double _getXmax(int, int) const;
  // double _getNumCoeffs(int) const;
  // double _getNPointsLeft(int) const;
  Eigen::MatrixXd& _getIntegralOperatorMatrix();
  Eigen::MatrixXd& _getInverseBornCSErrorMatrix();
  // Eigen::MatrixXd _permutation(int) const;
  // Eigen::MatrixXd _interpInvMatrix(int) const;
  // Eigen::MatrixXd _polConvKuraevFadinMatrix(int) const;
  // Eigen::MatrixXd _evalA(int) const;
  void _evalEqMatrix();
  // double _energySpreadWeight(double, std::size_t) const;
  Eigen::MatrixXd _energySpreadMatrix() const;
  // Eigen::RowVectorXd _polIntegralOp(int) const;
  // Eigen::RowVectorXd _splineIntegralOp() const;
  // Eigen::RowVectorXd _polGaussIntegralOp(int, int) const;
  // Eigen::RowVectorXd _evalPolGaussA(int, int) const;
  // Eigen::RowVectorXd _evalGaussDotProductOperator(std::size_t) const;
  // Eigen::RowVectorXd _evalPolA(int) const;
  const Eigen::RowVectorXd& _getDotProdOp() const;
  void _evalDotProductOperator();
  TF1* _createInterpFunction() const;
  TF1* _createDerivativeInterpFunction(std::size_t, const std::string&) const;
  // Eigen::RowVectorXd _interpProjector(double) const;
  // Eigen::RowVectorXd _interpDerivativeProjector(double, std::size_t) const;
  Eigen::VectorXd _bcsErr() const;
  // void _setDefaultInterpSettings();

 protected:
  Interpolator _interp;
 private:
  // std::vector<InterpPtSettings> _interpSettings;
  Eigen::MatrixXd _integralOperatorMatrix;
  Eigen::MatrixXd _invBornCSErrorMatrix;
  Eigen::RowVectorXd _dotProdOp;
};

#endif
