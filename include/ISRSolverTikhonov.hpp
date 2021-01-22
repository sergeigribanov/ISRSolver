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
  double evalEqNorm2() const;
  double evalEqNorm2NoErr() const;
  double evalApproxPerturbNorm2() const;
  double evalSmoothnessConstraintNorm2() const;
  double evalLCurveCurvature() const;
  double evalLCurveCurvatureDerivative() const;
  double evalApproxRegRelativeError(const Eigen::VectorXd&) const;
  double evalApproxPerturbRelativeError(const Eigen::VectorXd&,
					const Eigen::VectorXd&) const;
  void enableSolutionNorm2();
  void disableSolutionNorm2();
  void enableSolutionDerivativeNorm2();
  void disableSolutionDerivativeNorm2();

 protected:
  const Eigen::MatrixXd& _getInterpPointWiseDerivativeProjector() const;
  bool isSolutionNorm2Enabled() const;
  bool isSolutionNorm2DerivativeEnabled() const;
  double _evaldKsidAlpha(const Eigen::VectorXd&) const;
  double _evald2Ksid2Alpha(const Eigen::VectorXd&,
                           const Eigen::VectorXd&) const;
  void _evalProblemMatrices();
  void _evalInterpPointWiseDerivativeProjector();

 private:
  bool _enabledSolutionNorm2;
  bool _enabledSolutionDerivativeNorm2;
  double _alpha;
  Eigen::MatrixXd _interpPointWiseDerivativeProjector;
  Eigen::MatrixXd _mF;
  Eigen::MatrixXd _mR;
  Eigen::MatrixXd _mL;
  Eigen::FullPivLU<Eigen::MatrixXd> _luR;
  Eigen::FullPivLU<Eigen::MatrixXd> _luL;
};

#endif
