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

ISRSolverTikhonov::ISRSolverTikhonov(const std::string& inputPath,
                                     const InputOptions& inputOpts,
                                     double alpha)
    : ISRSolverSLAE(inputPath, inputOpts),
      _enabledSolutionNorm2(true),
      _enabledSolutionDerivativeNorm2(true),
      _alpha(alpha) {}

ISRSolverTikhonov::ISRSolverTikhonov(const ISRSolverTikhonov& solver) :
  ISRSolverSLAE::ISRSolverSLAE(solver),
  _enabledSolutionNorm2(solver._enabledSolutionNorm2),
  _enabledSolutionDerivativeNorm2(solver._enabledSolutionDerivativeNorm2),
  _alpha(solver._alpha),
  _interpPointWiseDerivativeProjector(solver._interpPointWiseDerivativeProjector)
{}

ISRSolverTikhonov::~ISRSolverTikhonov() {}

void ISRSolverTikhonov::solve() {
  if (!_isEqMatrixPrepared) {
    _evalDotProductOperator();
    _evalEqMatrix();
    _evalInterpPointWiseDerivativeProjector();
    _isEqMatrixPrepared = true;
  }
  _evalProblemMatrices();
  _bcs() = _luR.solve(_vcs());
  _getInverseBornCSErrorMatrix() = _mR.transpose() * _vcsInvErrMatrix() * _mR;
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

void ISRSolverTikhonov::_evalInterpPointWiseDerivativeProjector() {
  _interpPointWiseDerivativeProjector = Eigen::MatrixXd::Zero(_getN(), _getN());
  for (std::size_t i = 0; i < _getN(); ++i) {
    for (std::size_t j = 0; j < _getN(); ++j) {
      _interpPointWiseDerivativeProjector(i, j) = _interp.basisDerivEval(j, _ecm(i));
    }
  }
}

const Eigen::MatrixXd&
ISRSolverTikhonov::_getInterpPointWiseDerivativeProjector() const {
  return _interpPointWiseDerivativeProjector;
}

bool ISRSolverTikhonov::isSolutionNorm2Enabled() const {
  return _enabledSolutionNorm2;
}

bool ISRSolverTikhonov::isSolutionNorm2DerivativeEnabled() const {
  return _enabledSolutionDerivativeNorm2;
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
  return bcs().dot(_mF * bcs());
}

double ISRSolverTikhonov::evalLCurveCurvature() const {
  Eigen::VectorXd ds = -_luL.solve(bcs());
  double dksi = _evaldKsidAlpha(ds);
  return -std::fabs(1. / dksi / std::pow(1. + _alpha * _alpha, 1.5));
}

double ISRSolverTikhonov::evalLCurveCurvatureDerivative() const {
  Eigen::VectorXd ds = -_luL.solve(bcs());
  Eigen::VectorXd d2s = -2. * _luL.solve(ds);
  double dksi = _evaldKsidAlpha(ds);
  double d2ksi = _evald2Ksid2Alpha(ds, d2s);
  return -d2ksi * std::pow(dksi, -2.) * std::pow(1. + _alpha * _alpha, -1.5) +
      -3. * _alpha / dksi * std::pow(1. + _alpha * _alpha, -2.5);
}

double ISRSolverTikhonov::_evaldKsidAlpha(const Eigen::VectorXd& ds) const {
  return 2. * bcs().dot(_mF * ds);
}

double ISRSolverTikhonov::_evald2Ksid2Alpha(const Eigen::VectorXd& ds,
                                            const Eigen::VectorXd& d2s) const {
  return 2. * ds.dot(_mF * ds) + 2. * bcs().dot(_mF * d2s);
}

void ISRSolverTikhonov::_evalProblemMatrices() {
  Eigen::MatrixXd mAt = getIntegralOperatorMatrix().transpose();
  _mF = Eigen::MatrixXd::Zero(_getN(), _getN());
  if (isSolutionNorm2Enabled()) {
    _mF += _getDotProdOp().asDiagonal();
  }
  if (isSolutionNorm2DerivativeEnabled()) {
    _mF += _getInterpPointWiseDerivativeProjector().transpose() *
          (_getInterpPointWiseDerivativeProjector().array().colwise() *
           _getDotProdOp().transpose().array()).matrix();
  }
  _mR = getIntegralOperatorMatrix() +
        _alpha * _vcsInvErrMatrix().inverse() * mAt.inverse() * _mF;
  // std::cout << "---------------" << std::endl;
  // std::cout << "derivP:" << std::endl;
  // std::cout << _getInterpPointWiseDerivativeProjector() << std::endl;
  // std::cout << "_mF:" << std::endl;
  // std::cout << _mF << std::endl;
  // std::cout << "_mF.inverse():" << std::endl;
  // std::cout << _mF.inverse() << std::endl;
  // std::cout << "mAt:" << std::endl;
  // std::cout << mAt << std::endl;
  // std::cout << "_vcsInvErrMatrix():" << std::endl;
  // std::cout << _vcsInvErrMatrix() << std::endl;
  // std::cout << "getIntegralOperatorMatrix():" << std::endl;
  // std::cout << getIntegralOperatorMatrix() << std::endl;
  // std::cout << "_alpha = " << _alpha << std::endl;
  _mL = _mF.inverse() * mAt * _vcsInvErrMatrix() * getIntegralOperatorMatrix() +
        _alpha * Eigen::MatrixXd::Identity(_getN(), _getN());
  // std::cout << "_mL:" << std::endl;
  // std::cout << _mL << std::endl;
  _luR = Eigen::FullPivLU<Eigen::MatrixXd>(_mR);
  _luL = Eigen::FullPivLU<Eigen::MatrixXd>(_mL);
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
