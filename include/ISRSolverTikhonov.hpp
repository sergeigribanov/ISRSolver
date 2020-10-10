#ifndef __ISRSOLVER_TIKHONOV_HPP__
#define __ISRSOLVER_TIKHONOV_HPP__

#include "ISRSolverSLAE.hpp"

class ISRSolverTikhonov : public ISRSolverSLAE {
 public:
  ISRSolverTikhonov(const std::string& inputPath, const InputOptions& inputOpts,
                    double alpha = 1.e-2);
  virtual ~ISRSolverTikhonov();
  virtual void solve() override;
  double getAlpha() const;
  const Eigen::MatrixXd& getHessian() const;
  void enableSolutionNorm2();
  void disableSolutionNorm2();
  void enableSolutionDerivativeNorm2();
  void disableSolutionDerivativeNorm2();
  void setAlpha(double);

 protected:
  Eigen::RowVectorXd _polIntegralOp(int) const;
  Eigen::RowVectorXd _evalPolA(int) const;
  const Eigen::RowVectorXd& _getDotProdOp() const;
  const Eigen::MatrixXd& _getInterpPointWiseDerivativeProjector() const;
  bool isSolutionNorm2Enabled() const;
  bool isSolutionNorm2DerivativeEnabled() const;
  double _evalRegFuncNorm2(const Eigen::VectorXd&) const;
  Eigen::VectorXd _evalRegFuncGradNorm2(const Eigen::VectorXd&) const;
  void _evalHessian();
  void _evalDotProductOperator();
  void _evalInterpPointWiseDerivativeProjector();

 private:
  bool _enabledSolutionNorm2;
  bool _enabledSolutionDerivativeNorm2;
  double _alpha;
  Eigen::RowVectorXd _dotProdOp;
  Eigen::MatrixXd _interpPointWiseDerivativeProjector;
  Eigen::MatrixXd _hessian;
  friend double objectiveFCN(unsigned, const double*, double*, void*);
};

#endif
