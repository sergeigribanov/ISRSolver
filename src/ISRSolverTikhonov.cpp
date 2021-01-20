#include "ISRSolverTikhonov.hpp"

#include <TFile.h>
#include <TGraphErrors.h>
#include <TMatrixD.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <nlohmann/json.hpp>
#include <nlopt.hpp>
#include <set>

#include "integration.hpp"

double objectiveFCN(unsigned n, const double* z, double* grad, void* solver) {
  auto sp = reinterpret_cast<ISRSolverTikhonov*>(solver);
  Eigen::Map<const Eigen::VectorXd> vm(z, n);
  if (grad) {
    Eigen::VectorXd vgrad = sp->_evalRegFuncGradNorm2(vm);
    for (std::size_t i = 0; i < n; ++i) {
      grad[i] = vgrad(i);
    }
  }
  return sp->_evalRegFuncNorm2(vm);
}

ISRSolverTikhonov::ISRSolverTikhonov(const std::string& inputPath,
                                     const InputOptions& inputOpts,
                                     double alpha)
    : ISRSolverSLAE(inputPath, inputOpts),
      _solutionPositivity(false),
      _enabledSolutionNorm2(true),
      _enabledSolutionDerivativeNorm2(true),
      _alpha(alpha) {}

ISRSolverTikhonov::ISRSolverTikhonov(const ISRSolverTikhonov& solver) :
  ISRSolverSLAE::ISRSolverSLAE(solver),
  _solutionPositivity(solver._solutionPositivity),
  _enabledSolutionNorm2(solver._enabledSolutionNorm2),
  _enabledSolutionDerivativeNorm2(solver._enabledSolutionDerivativeNorm2),
  _alpha(solver._alpha),
  _interpPointWiseDerivativeProjector(solver._interpPointWiseDerivativeProjector),
  _hessian(solver._hessian) {}

ISRSolverTikhonov::~ISRSolverTikhonov() {}

void ISRSolverTikhonov::solve() {
  _evalDotProductOperator();
  _evalEqMatrix();
  _evalInterpPointWiseDerivativeProjector();
  _evalHessian();
  nlopt::opt opt(nlopt::LD_MMA, _getN());
  if (_solutionPositivity) {
    std::vector<double> lowerBounds(_getN(), 0);
    opt.set_lower_bounds(lowerBounds);
  }
  opt.set_min_objective(objectiveFCN, this);
  opt.set_xtol_rel(1.e-6);
  std::vector<double> z(_getN());
  std::transform(_vcs().data(), _vcs().data() + _getN(), z.begin(),
		 [](double arg) {return std::max(0., arg);});
  double minf;
  try {
    opt.optimize(z, minf);
    _bcs() = Eigen::Map<Eigen::VectorXd>(z.data(), _getN());
    _getInverseBornCSErrorMatrix() = 0.5 * getHessian();
  } catch (std::exception& e) {
    std::cerr << "nlopt failed: " << e.what() << std::endl;
  }
}

double ISRSolverTikhonov::getAlpha() const { return _alpha; }

void ISRSolverTikhonov::enableSolutionNorm2() { _enabledSolutionNorm2 = true; }

void ISRSolverTikhonov::disableSolutionNorm2() {
  _enabledSolutionNorm2 = false;
}

void ISRSolverTikhonov::enableSolutionDerivativeNorm2() {
  _enabledSolutionDerivativeNorm2 = true;
}

void ISRSolverTikhonov::disableSolutionDerivativeNorm2() {
  _enabledSolutionDerivativeNorm2 = false;
}

void ISRSolverTikhonov::setAlpha(double alpha) { _alpha = alpha; }

double ISRSolverTikhonov::_evalRegFuncNorm2(const Eigen::VectorXd& z) const {
  Eigen::VectorXd dv = getIntegralOperatorMatrix() * z - _vcs();
  Eigen::VectorXd dz = _getInterpPointWiseDerivativeProjector() * z;
  double result = dv.dot(_vcsInvErrMatrix() * dv);
  if (isSolutionNorm2Enabled()) {
    result += _alpha * _getDotProdOp() * (z.array() * z.array()).matrix();
  }
  if (isSolutionNorm2DerivativeEnabled()) {
    result += _alpha * _getDotProdOp() * (dz.array() * dz.array()).matrix();
  }
  return result;
}

void ISRSolverTikhonov::_evalInterpPointWiseDerivativeProjector() {
  _interpPointWiseDerivativeProjector = Eigen::MatrixXd::Zero(_getN(), _getN());
  for (std::size_t j = 0; j < _getN(); ++j) {
    Eigen::MatrixXd coeffs = _interpInvMatrix(j) * _permutation(j);
    const std::size_t nc = _getNumCoeffs(j);
    for (std::size_t k = 1; k < nc; ++k) {
      _interpPointWiseDerivativeProjector.row(j) +=
          coeffs.row(k) * k * std::pow(_s(j), k - 1);
    }
  }
}

const Eigen::MatrixXd&
ISRSolverTikhonov::_getInterpPointWiseDerivativeProjector() const {
  return _interpPointWiseDerivativeProjector;
}

Eigen::VectorXd ISRSolverTikhonov::_evalRegFuncGradNorm2(
    const Eigen::VectorXd& z) const {
  Eigen::VectorXd dv = getIntegralOperatorMatrix() * z - _vcs();
  Eigen::VectorXd dz = _getInterpPointWiseDerivativeProjector() * z;
  Eigen::VectorXd result = 2. * getIntegralOperatorMatrix().transpose() *
                           _vcsInvErrMatrix() * dv;
  if (isSolutionNorm2Enabled()) {
    result += 2. * _alpha * (_getDotProdOp().transpose().array() * z.array()).matrix();
  }
  if (isSolutionNorm2DerivativeEnabled()) {
    result += 2. * _alpha *
              _getInterpPointWiseDerivativeProjector().transpose() *
              (_getDotProdOp().transpose().array() * dz.array()).matrix();
  }
  return result;
}

bool ISRSolverTikhonov::isSolutionNorm2Enabled() const {
  return _enabledSolutionNorm2;
}

bool ISRSolverTikhonov::isSolutionNorm2DerivativeEnabled() const {
  return _enabledSolutionDerivativeNorm2;
}

void ISRSolverTikhonov::_evalHessian() {
  _hessian = 2. * getIntegralOperatorMatrix().transpose() *
             _vcsInvErrMatrix() *
             getIntegralOperatorMatrix();
  if (isSolutionNorm2Enabled()) {
    _hessian += 2. * _alpha * _getDotProdOp().asDiagonal();
  }
  if (isSolutionNorm2DerivativeEnabled()) {
   _hessian += 2. * _alpha *
               _getInterpPointWiseDerivativeProjector().transpose() *
               (_getInterpPointWiseDerivativeProjector().array().colwise() *
                _getDotProdOp().transpose().array()).matrix();
  }
}

const Eigen::MatrixXd& ISRSolverTikhonov::getHessian() const {
  return _hessian;
}

void ISRSolverTikhonov::disableSolutionPositivity() {
  _solutionPositivity = false;
}

void ISRSolverTikhonov::enableSolutionPositivity() {
  _solutionPositivity = true;
}

double ISRSolverTikhonov::evalEqNorm2() const {
  Eigen::VectorXd dv = getIntegralOperatorMatrix() * bcs() - _vcs();
  return dv.dot(_vcsInvErrMatrix() * dv);
}

double ISRSolverTikhonov::evalEqNorm2NoErr() const {
  Eigen::VectorXd dv = getIntegralOperatorMatrix() * bcs() - _vcs();
  return dv.dot(dv);
}

double ISRSolverTikhonov::evalApproxPerturbNorm2() const {
  return _getDotProdOp() * _vcsErr().array().pow(2.).matrix();
}

double ISRSolverTikhonov::evalSmoothnessConstraintNorm2() const {
  Eigen::VectorXd dz = _getInterpPointWiseDerivativeProjector() * bcs();
  double result = 0;
  if (isSolutionNorm2Enabled()) {
    result += _getDotProdOp() * bcs().array().pow(2.).matrix();
  }
  if (isSolutionNorm2DerivativeEnabled()) {
    result += _getDotProdOp() * dz.array().pow(2.).matrix();
  }
  return result;
}

double ISRSolverTikhonov::evalLCurveCurvature() const {
  Eigen::MatrixXd mQ;
  Eigen::MatrixXd mF;
  std::tie(mQ, mF) = _evaldKsiMatrices();
  Eigen::VectorXd ds = -mQ * mF * mQ *
                       getIntegralOperatorMatrix().transpose() *
                       _vcsInvErrMatrix() *
                       _vcs();
  double dksi = _evaldKsidAlpha(ds);
  return -std::fabs(1. / dksi / std::pow(1. + _alpha * _alpha, 1.5));
}

double ISRSolverTikhonov::evalUCurve() const {
  double rho = evalEqNorm2();
  double ksi = evalSmoothnessConstraintNorm2();
  return 1. / rho + 1. / ksi;
}

double ISRSolverTikhonov::evalUCurveCurvature() const {
  Eigen::MatrixXd mQ;
  Eigen::MatrixXd mF;
  std::tie(mQ, mF) = _evaldKsiMatrices();
  Eigen::VectorXd ds = -mQ * mF * mQ *
                       getIntegralOperatorMatrix().transpose() *
                       _vcsInvErrMatrix() *
                       _vcs();
  Eigen::VectorXd d2s = -2. * mQ * mF * ds;
  double dksi = _evaldKsidAlpha(ds);
  double d2ksi = _evald2Ksid2Alpha(ds, d2s);
  double rho = evalEqNorm2();
  double ksi = evalSmoothnessConstraintNorm2();
  double du = dksi * (_alpha / std::pow(rho, 2) - 1. / std::pow(ksi, 2));
  double d2u = d2ksi * (_alpha / std::pow(rho, 2) - 1. / std::pow(ksi, 2)) +
               dksi * (1. / std::pow(rho, 2) +
                       2. * std::pow(_alpha, 2) * dksi / std::pow(rho, 3) +
                       2. * dksi / std::pow(ksi, 3));
  return -std::fabs(d2u / std::sqrt(1. + du * du));
}

double ISRSolverTikhonov::_evaldKsidAlpha(const Eigen::VectorXd& ds) const {
  Eigen::VectorXd dz = _getInterpPointWiseDerivativeProjector() * bcs();
  Eigen::VectorXd dsz = _getInterpPointWiseDerivativeProjector() * ds;
  double result = 0;
  if (isSolutionNorm2Enabled()) {
    result += 2. * _getDotProdOp() * (bcs().array() * ds.array()).matrix();
  }
  if (isSolutionNorm2DerivativeEnabled()) {
    result += 2. * _getDotProdOp() * (dz.array() * dsz.array()).matrix();
  }
  return result;
}

double ISRSolverTikhonov::_evald2Ksid2Alpha(const Eigen::VectorXd& ds,
                                          const Eigen::VectorXd& d2s) const {
  Eigen::VectorXd dz = _getInterpPointWiseDerivativeProjector() * bcs();
  Eigen::VectorXd dsz = _getInterpPointWiseDerivativeProjector() * ds;
  Eigen::VectorXd d2sz = _getInterpPointWiseDerivativeProjector() * d2s;
  double result = 0;
  if (isSolutionNorm2Enabled()) {
    result += 2. * _getDotProdOp() * (ds.array().pow(2.) +
                                      bcs().array() * d2s.array()).matrix();
  }
  if (isSolutionNorm2DerivativeEnabled()) {
    result += _getDotProdOp() * (dsz.array().pow(2.) +
                                 dz.array() * d2sz.array()).matrix();
  }
  return result;
}

std::pair<Eigen::MatrixXd, Eigen::MatrixXd>
ISRSolverTikhonov::_evaldKsiMatrices() const {
  Eigen::MatrixXd mF = Eigen::MatrixXd::Zero(_getN(), _getN());
  Eigen::MatrixXd mAt = getIntegralOperatorMatrix().transpose();
  if (isSolutionNorm2Enabled()) {
    mF += _getDotProdOp().asDiagonal();
  }
  if (isSolutionNorm2DerivativeEnabled()) {
    mF += _getInterpPointWiseDerivativeProjector().transpose() *
          (_getInterpPointWiseDerivativeProjector().array().colwise() *
           _getDotProdOp().transpose().array()).matrix();
  }
  Eigen::MatrixXd mQ = mAt * _vcsInvErrMatrix() * getIntegralOperatorMatrix() +
                       _alpha * mF;
  mQ = mQ.inverse();
  return std::make_pair(mQ, mF);
}

// !!! TO DO: modify reg and perturb errors
double ISRSolverTikhonov::evalApproxRegRelativeError(const Eigen::VectorXd& origBCS) const {
  ISRSolverTikhonov solver(*this);
  solver._vcs() = getIntegralOperatorMatrix() * origBCS;
  solver.solve();
  Eigen::VectorXd dbcs = (origBCS - solver.bcs());
  double errNorm = _getDotProdOp() * (dbcs.array() * dbcs.array()).matrix();
  errNorm = std::sqrt(errNorm);
  double bcsNorm = _getDotProdOp() * (origBCS.array() * origBCS.array()).matrix();
  bcsNorm = std::sqrt(bcsNorm);
  return errNorm / bcsNorm;
}

double ISRSolverTikhonov::evalApproxPerturbRelativeError(const Eigen::VectorXd& origBCS,
							 const Eigen::VectorXd& vcsPerturbation) const {
  ISRSolverTikhonov solver(*this);
  solver._vcs() = vcsPerturbation;
  solver.solve();
  double errNorm = _getDotProdOp() * (solver.bcs().array() * solver.bcs().array()).matrix();
  errNorm = std::sqrt(errNorm);
  double bcsNorm = _getDotProdOp() * (origBCS.array() * origBCS.array()).matrix();
  bcsNorm = std::sqrt(bcsNorm);
  return errNorm / bcsNorm;
}
