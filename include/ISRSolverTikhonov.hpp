#ifndef __ISRSOLVER_TIKHONOV_HPP__
#define __ISRSOLVER_TIKHONOV_HPP__

#include "ISRSolverSLAE.hpp"

class ISRSolverTikhonov : public ISRSolverSLAE {
 public:
  ISRSolverTikhonov(const std::string& inputPath, const InputOptions& inputOpts,
                    double alpha = 1.e-2);
  ISRSolverTikhonov(const ISRSolverTikhonov&);
  virtual ~ISRSolverTikhonov();
  virtual void solve() override;
  double getAlpha() const;
  void setAlpha(double);
  const Eigen::MatrixXd& getHessian() const;
  double evalEqNorm2() const;
  double evalSmoothnessConstraintNorm2() const;
  double evalCurvature() const;
  double evalApproxRegRelativeError(const Eigen::VectorXd&) const;
  double evalApproxPerturbRelativeError(const Eigen::VectorXd&,
					const Eigen::VectorXd&) const;
  void enableSolutionNorm2();
  void disableSolutionNorm2();
  void enableSolutionDerivativeNorm2();
  void disableSolutionDerivativeNorm2();
  void disableSolutionPositivity();
  void enableSolutionPositivity();

 protected:
  Eigen::RowVectorXd _polIntegralOp(int) const;
  Eigen::RowVectorXd _evalPolA(int) const;
  const Eigen::RowVectorXd& _getDotProdOp() const;
  const Eigen::MatrixXd& _getInterpPointWiseDerivativeProjector() const;
  bool isSolutionNorm2Enabled() const;
  bool isSolutionNorm2DerivativeEnabled() const;
  double _evalRegFuncNorm2(const Eigen::VectorXd&) const;
  Eigen::VectorXd _evalRegFuncGradNorm2(const Eigen::VectorXd&) const;
  double _evaldKsidAlpha() const;
  Eigen::VectorXd _evaldSoldAlpha() const;
  void _evalHessian();
  void _evalDotProductOperator();
  void _evalInterpPointWiseDerivativeProjector();

 private:
  bool _solutionPositivity;
  bool _enabledSolutionNorm2;
  bool _enabledSolutionDerivativeNorm2;
  double _alpha;
  Eigen::RowVectorXd _dotProdOp;
  Eigen::MatrixXd _interpPointWiseDerivativeProjector;
  Eigen::MatrixXd _hessian;
  friend double objectiveFCN(unsigned, const double*, double*, void*);
};

#endif
