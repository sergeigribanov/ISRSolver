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
  _dotProdOp(solver._dotProdOp),
  _interpPointWiseDerivativeProjector(solver._interpPointWiseDerivativeProjector),
  _hessian(solver._hessian) {}

ISRSolverTikhonov::~ISRSolverTikhonov() {}

void ISRSolverTikhonov::solve() {
  _evalEqMatrix();
  _evalDotProductOperator();
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
  Eigen::VectorXd invVCSErr2 = _vcsErr().array().pow(-2.).matrix();
  Eigen::VectorXd dz = _getInterpPointWiseDerivativeProjector() * z;
  return _getDotProdOp() *
         (dv.array() * dv.array() * invVCSErr2.array() +
          _enabledSolutionNorm2 * _alpha * z.array() * z.array() +
          _enabledSolutionDerivativeNorm2 * _alpha * dz.array() * dz.array())
             .matrix();
}

Eigen::RowVectorXd ISRSolverTikhonov::_polIntegralOp(int j) const {
  const std::size_t nc = _getNumCoeffs(j);
  Eigen::VectorXd result(nc);
  std::size_t k;
  double error;
  double sMin;
  if (j == 0) {
    sMin = _sThreshold();
  } else {
    sMin = _s(j - 1);
  }
  double sMax = _s(j);
  std::function<double(double)> fcn = [&k](double t) { return std::pow(t, k); };
  for (k = 0; k < nc; ++k) {
    result(k) = integrate(fcn, sMin, sMax, error);
  }
  return result.transpose();
}

Eigen::RowVectorXd ISRSolverTikhonov::_evalPolA(int j) const {
  return _polIntegralOp(j) * _interpInvMatrix(j) * _permutation(j);
}

void ISRSolverTikhonov::_evalDotProductOperator() {
  _dotProdOp = Eigen::RowVectorXd::Zero(_getN());
  for (std::size_t i = 0; i < _getN(); ++i) {
    _dotProdOp += _evalPolA(i);
  }
}

const Eigen::RowVectorXd& ISRSolverTikhonov::_getDotProdOp() const {
  return _dotProdOp;
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
  Eigen::VectorXd invVCSErr2 = _vcsErr().array().pow(-2.).matrix();
  Eigen::VectorXd result = Eigen::VectorXd::Zero(_getN());
  Eigen::VectorXd dz = _getInterpPointWiseDerivativeProjector() * z;
  for (std::size_t l = 0; l < _getN(); ++l) {
    for (std::size_t i = 0; i < _getN(); ++i) {
      result(l) += 2 * _getDotProdOp()(i) * getIntegralOperatorMatrix()(i, l) *
                   dv(i) * invVCSErr2(i);
      if (isSolutionNorm2DerivativeEnabled()) {
        result(l) += 2 * _alpha * _getDotProdOp()(i) *
                     _getInterpPointWiseDerivativeProjector()(i, l) * dz(i);
      }
      if (l == i && isSolutionNorm2Enabled()) {
        result(l) += 2 * _alpha * _getDotProdOp()(i) * z(i);
      }
    }
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
  Eigen::VectorXd invVCSErr2 = _vcsErr().array().pow(-2.).matrix();
  _hessian = Eigen::MatrixXd::Zero(_getN(), _getN());
  for (std::size_t l = 0; l < _getN(); ++l) {
    for (std::size_t m = 0; m < _getN(); ++m) {
      for (std::size_t i = 0; i < _getN(); ++i) {
        _hessian(l, m) += 2 * getIntegralOperatorMatrix()(i, l) *
                          getIntegralOperatorMatrix()(i, m) * invVCSErr2(i);
        if (isSolutionNorm2DerivativeEnabled()) {
          _hessian(l, m) += 2 * _alpha *
                            _getInterpPointWiseDerivativeProjector()(i, l) *
                            _getInterpPointWiseDerivativeProjector()(i, m);
        }
        if (l == i && m == i && isSolutionNorm2Enabled()) {
          _hessian(l, m) += 2 * _alpha;
        }
      }
    }
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
  Eigen::VectorXd invVCSErr2 = _vcsErr().array().pow(-2.).matrix();
  return _getDotProdOp() *
    (dv.array() * dv.array() * invVCSErr2.array()).matrix();
}

double ISRSolverTikhonov::evalEqNorm2NoErr() const {
  Eigen::VectorXd dv = getIntegralOperatorMatrix() * bcs() - _vcs();
  return _getDotProdOp() *
    (dv.array() * dv.array()).matrix();
}

double ISRSolverTikhonov::evalApproxPerturbNorm2() const {
  return _getDotProdOp() *
    (_vcsErr().array() * _vcsErr().array()).matrix();
}

double ISRSolverTikhonov::evalSmoothnessConstraintNorm2() const {
  Eigen::VectorXd dv = getIntegralOperatorMatrix() * bcs() - _vcs();
  Eigen::VectorXd invVCSErr2 = _vcsErr().array().pow(-2.).matrix();
  Eigen::VectorXd dz = _getInterpPointWiseDerivativeProjector() * bcs();
  return _getDotProdOp() *
    (_enabledSolutionNorm2 * bcs().array() * bcs().array() +
     _enabledSolutionDerivativeNorm2 * dz.array() * dz.array()).matrix();
}

double ISRSolverTikhonov::evalCurvature() const {
  double ksi = evalSmoothnessConstraintNorm2();
  double rho = evalEqNorm2();
  double dksi = _evaldKsidAlpha();
  return 2. * ksi * rho / dksi * (rho * ksi + _alpha * _alpha * ksi * dksi + _alpha * rho * dksi) /
    std::pow(_alpha * _alpha * ksi * ksi + rho * rho, 1.5);
}

double ISRSolverTikhonov::_evaldKsidAlpha() const {
  Eigen::VectorXd ds = _evaldSoldAlpha();
  Eigen::VectorXd dz = _getInterpPointWiseDerivativeProjector() * bcs();
  Eigen::VectorXd dsz = _getInterpPointWiseDerivativeProjector() * ds;
  return _getDotProdOp() * (_enabledSolutionNorm2 * bcs().array() * ds.array() +
			    _enabledSolutionDerivativeNorm2 * dz.array() * dsz.array()).matrix();
}

Eigen::VectorXd ISRSolverTikhonov::_evaldSoldAlpha() const {
  Eigen::MatrixXd mQ = Eigen::MatrixXd::Zero(_getN(), _getN());
  Eigen::MatrixXd mF = Eigen::MatrixXd::Zero(_getN(), _getN());
  Eigen::MatrixXd mL = Eigen::MatrixXd::Zero(_getN(), _getN());
  Eigen::VectorXd vB = Eigen::VectorXd::Zero(_getN());
  Eigen::VectorXd invVCSErr2 = _vcsErr().array().pow(-2.).matrix();
  for (std::size_t i = 0; i < _getN(); ++i) {
    for (std::size_t l = 0; l < _getN(); ++l) {
      vB(l) += _getDotProdOp()(i) * getIntegralOperatorMatrix()(i, l) * _vcs()(i) * invVCSErr2(i);
      if (i == l && _enabledSolutionNorm2) {
	mF(l, i) += _getDotProdOp()(i) * bcs()(i);
      }
      for (std::size_t k = 0; k < _getN(); ++k) {
	mQ(l, k) += _getDotProdOp()(i) * getIntegralOperatorMatrix()(i, l) *
	  getIntegralOperatorMatrix()(i, k) * invVCSErr2(i);
	if (_enabledSolutionDerivativeNorm2) {
	  mL(l, k) += _getDotProdOp()(i) * _getInterpPointWiseDerivativeProjector()(i, l) *
	    _getInterpPointWiseDerivativeProjector()(i, k);
	}
      }
    }
  }
  Eigen::MatrixXd mP = mF + mL;
  Eigen::MatrixXd mT = mQ + _alpha * mP;
  mT = mT.inverse();
  return -mT * mP * mT * vB;
}

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
