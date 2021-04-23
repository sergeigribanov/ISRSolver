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

#include "Integration.hpp"

ISRSolverTikhonov::ISRSolverTikhonov(TGraphErrors* vcsGraph,
                                     double thresholdEnergy,
                                     double lambda) :
    ISRSolverSLE(vcsGraph, thresholdEnergy),
    _enabledDerivNorm2Reg(true),
    _lambda(lambda) {}

ISRSolverTikhonov::ISRSolverTikhonov(TGraphErrors* vcsGraph,
                                     TEfficiency* eff,
                                     double thresholdEnergy,
                                     double lambda) :
    ISRSolverSLE(vcsGraph, eff, thresholdEnergy),
    _enabledDerivNorm2Reg(true),
    _lambda(lambda) {}

ISRSolverTikhonov::ISRSolverTikhonov(const std::string& inputPath,
                                     const InputOptions& inputOpts,
                                     double lambda)
    : ISRSolverSLE(inputPath, inputOpts),
      _enabledDerivNorm2Reg(true),
      _lambda(lambda) {}

ISRSolverTikhonov::ISRSolverTikhonov(const ISRSolverTikhonov& solver) :
    ISRSolverSLE(solver),
    _enabledDerivNorm2Reg(true),
    _lambda(solver._lambda),
    _interpPointWiseDerivativeProjector(solver._interpPointWiseDerivativeProjector) {}

ISRSolverTikhonov::~ISRSolverTikhonov() {}

void ISRSolverTikhonov::solve() {
  if (!_isEqMatrixPrepared) {
    _evalDotProductOperator();
    evalEqMatrix();
    _evalInterpPointWiseDerivativeProjector();
    _isEqMatrixPrepared = true;
  }
  _evalProblemMatrices();
  Eigen::MatrixXd mAp = _luT.solve(getIntegralOperatorMatrix().transpose() * _vcsInvCovMatrix());
  _bcs() = mAp * _vcs();
  _getBornCSCovMatrix() = mAp * _vcsInvCovMatrix().inverse() * mAp.transpose();
}

double ISRSolverTikhonov::getLambda() const {
  return _lambda;
}

bool ISRSolverTikhonov::isDerivNorm2RegIsEnabled() const {
  return _enabledDerivNorm2Reg;
}

void ISRSolverTikhonov::enableDerivNorm2Regularizator() {
  _enabledDerivNorm2Reg = true;
}

void ISRSolverTikhonov::disableDerivNorm2Regularizator() {
  _enabledDerivNorm2Reg = false;
}

void ISRSolverTikhonov::setLambda(double lambda) { _lambda = lambda; }

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

double ISRSolverTikhonov::evalEqNorm2() const {
  Eigen::VectorXd dv = getIntegralOperatorMatrix() * bcs() - _vcs();
  return dv.dot(_vcsInvCovMatrix() * dv);
}

double ISRSolverTikhonov::evalSmoothnessConstraintNorm2() const {
  return bcs().dot(_mF * bcs());
}

double ISRSolverTikhonov::evalLCurveCurvature() const {
  Eigen::VectorXd ds = -_luL.solve(bcs());
  double dksi = _evaldKsidLambda(ds);
  return -std::fabs(1. / dksi / std::pow(1. + _lambda * _lambda, 1.5));
}

double ISRSolverTikhonov::evalLCurveCurvatureDerivative() const {
  Eigen::VectorXd ds = -_luL.solve(bcs());
  Eigen::VectorXd d2s = -2. * _luL.solve(ds);
  double dksi = _evaldKsidLambda(ds);
  double d2ksi = _evald2Ksid2Lambda(ds, d2s);
  return -d2ksi * std::pow(dksi, -2.) * std::pow(1. + _lambda * _lambda, -1.5) +
      -3. * _lambda / dksi * std::pow(1. + _lambda * _lambda, -2.5);
}

double ISRSolverTikhonov::_evaldKsidLambda(const Eigen::VectorXd& ds) const {
  return 2. * bcs().dot(_mF * ds);
}

double ISRSolverTikhonov::_evald2Ksid2Lambda(const Eigen::VectorXd& ds,
                                            const Eigen::VectorXd& d2s) const {
  return 2. * ds.dot(_mF * ds) + 2. * bcs().dot(_mF * d2s);
}

void ISRSolverTikhonov::_evalProblemMatrices() {
  Eigen::MatrixXd mAt = getIntegralOperatorMatrix().transpose();
  _mF = Eigen::MatrixXd::Zero(_getN(), _getN());
  if (isDerivNorm2RegIsEnabled()) {
    _mF += _getInterpPointWiseDerivativeProjector().transpose() *
           (_getInterpPointWiseDerivativeProjector().array().colwise() *
            _getDotProdOp().transpose().array()).matrix();
  } else {
    _mF += _getDotProdOp().asDiagonal();
  }
  Eigen::MatrixXd mT = mAt * _vcsInvCovMatrix() * getIntegralOperatorMatrix() +
        _lambda * _mF;
  _mL = _mF.inverse() * mAt * _vcsInvCovMatrix() * getIntegralOperatorMatrix() +
        _lambda * Eigen::MatrixXd::Identity(_getN(), _getN());
  _luT = Eigen::FullPivLU<Eigen::MatrixXd>(mT);
  _luL = Eigen::FullPivLU<Eigen::MatrixXd>(_mL);
}
