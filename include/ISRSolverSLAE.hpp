#ifndef __ISRSOLVER_SLAE_HPP__
#define __ISRSOLVER_SLAE_HPP__

#include <TF1.h>

#include "BaseISRSolver.hpp"

class ISRSolverSLAE : public BaseISRSolver {
 public:
  ISRSolverSLAE(const std::string& inputPath, const InputOptions& inputOpts);
  virtual ~ISRSolverSLAE();
  const Eigen::MatrixXd& getIntegralOperatorMatrix() const;
  const Eigen::MatrixXd& getInverseBornCSErrorMatrix() const;
  virtual void solve() override;
  virtual void save(const std::string& outputPath,
                    const OutputOptions& outputOpts) override;
  void setInterpSettings(const std::vector<InterpPtSettings>&) noexcept(false);
  void setInterpSettings(const std::string&) noexcept(false);

 protected:
  double _getXmin(int, int) const;
  double _getXmax(int, int) const;
  double _getNumCoeffs(int) const;
  double _getNPointsLeft(int) const;
  Eigen::MatrixXd& _getIntegralOperatorMatrix();
  Eigen::MatrixXd& _getInverseBornCSErrorMatrix();
  Eigen::MatrixXd _permutation(int) const;
  Eigen::MatrixXd _interpInvMatrix(int) const;
  Eigen::MatrixXd _polConvKuraevFadinMatrix(int) const;
  Eigen::MatrixXd _evalA(int) const;
  void _evalEqMatrix();
  TF1* _createInterpFunction() const;
  TF1* _createDerivativeInterpFunction(std::size_t, const std::string&) const;
  Eigen::RowVectorXd _interpProjector(double) const;
  Eigen::RowVectorXd _interpDerivativeProjector(double, std::size_t) const;
  Eigen::VectorXd _bcsErr() const;
  void _setDefaultInterpSettings();

 private:
  std::vector<InterpPtSettings> _interpSettings;
  Eigen::MatrixXd _integralOperatorMatrix;
  Eigen::MatrixXd _invBornCSErrorMatrix;
};

#endif
