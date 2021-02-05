#include "ISRSolverSLAE.hpp"

#include <TFile.h>
#include <TGraphErrors.h>
#include <TMatrixD.h>

#include <algorithm>
#include <boost/lexical_cast.hpp>
#include <cmath>
#include <fstream>
#include <nlohmann/json.hpp>
#include <set>

#include <iostream>

#include "integration.hpp"
#include "kuraev_fadin.hpp"

using json = nlohmann::json;

ISRSolverSLAE::ISRSolverSLAE(const std::string& inputPath,
                             const InputOptions& inputOpts) :
    BaseISRSolver(inputPath, inputOpts),
    _interp(Interpolator(ecm(), getThresholdEnergy())),
    _isEqMatrixPrepared(false) {}

ISRSolverSLAE::ISRSolverSLAE(const ISRSolverSLAE& solver) :
  BaseISRSolver::BaseISRSolver(solver),
  _interp(solver._interp),
  _isEqMatrixPrepared(solver._isEqMatrixPrepared),
  _integralOperatorMatrix(solver._integralOperatorMatrix),
  _covMatrixBornCS(solver._covMatrixBornCS),
  _dotProdOp(solver._dotProdOp) {}

ISRSolverSLAE::~ISRSolverSLAE() {}

const Eigen::MatrixXd& ISRSolverSLAE::getIntegralOperatorMatrix() const {
  return _integralOperatorMatrix;
}

const Eigen::MatrixXd& ISRSolverSLAE::getBornCSCovMatrix() const {
  return _covMatrixBornCS;
}

Eigen::MatrixXd& ISRSolverSLAE::_getIntegralOperatorMatrix() {
  return _integralOperatorMatrix;
}

Eigen::MatrixXd& ISRSolverSLAE::_getBornCSCovMatrix() {
  return _covMatrixBornCS;
}

void ISRSolverSLAE::solve() {
  if (!_isEqMatrixPrepared) {
    _evalEqMatrix();
    _isEqMatrixPrepared = true;
  }
  _bcs() =
      _integralOperatorMatrix.completeOrthogonalDecomposition().solve(_vcs());
  _covMatrixBornCS = (_integralOperatorMatrix.transpose() *
                      _vcsInvCovMatrix() * _integralOperatorMatrix).inverse();
}

void ISRSolverSLAE::save(const std::string& outputPath,
                         const OutputOptions& outputOpts) {
  TGraphErrors vcs(_getN(), _ecm().data(), _vcs().data(), _ecmErr().data(),
                   _vcsErr().data());
  TGraphErrors bcs(_getN(), _ecm().data(), _bcs().data(), 0, _bcsErr().data());
  TMatrixD intergalOperatorMatrix(_getN(), _getN());
  Eigen::MatrixXd tmpIntOpM = _integralOperatorMatrix.transpose();
  intergalOperatorMatrix.SetMatrixArray(tmpIntOpM.data());
  TMatrixD bornCSCovMatrix(_getN(), _getN());
  Eigen::MatrixXd tmpCovM = _covMatrixBornCS.transpose();
  bornCSCovMatrix.SetMatrixArray(tmpCovM.data());
  auto f0 = _createInterpFunction();
  auto f1 = _createDerivativeInterpFunction(1, "interp1DivFCN");
  auto fl = TFile::Open(outputPath.c_str(), "recreate");
  fl->cd();
  vcs.Write(outputOpts.visibleCSGraphName.c_str());
  bcs.Write(outputOpts.bornCSGraphName.c_str());
  intergalOperatorMatrix.Write("intergalOperatorMatrix");
  bornCSCovMatrix.Write("covMatrixBornCS");
  f0->Write();
  f1->Write();
  fl->Close();
  delete f0;
  delete f1;
  delete fl;
}

void ISRSolverSLAE::setRangeInterpSettings(
    const std::vector<std::tuple<bool, int, int>>& interpRangeSettings) {
  _interp = Interpolator(interpRangeSettings, ecm(), getThresholdEnergy());
}

void ISRSolverSLAE::setRangeInterpSettings(const std::string& pathToJSON) {
  _interp = Interpolator(pathToJSON, ecm(), getThresholdEnergy());
}

void ISRSolverSLAE::_evalEqMatrix() {
  std::size_t j;
  std::function<double(double)> fcn =
      [&j, this](double energy) {
        double result = this->_interp.basisEval(j, energy);
        return result;
      };
  _integralOperatorMatrix = Eigen::MatrixXd::Zero(_getN(), _getN());
  // TO DO: optimize
  for (j = 0; j < _getN(); ++j) {
    for(std::size_t i = 0; i < _getN(); ++i) {
      const double x_max =
        1. - std::pow(getThresholdEnergy() / _ecm(i), 2);
      _integralOperatorMatrix(i, j) =
          kuraev_fadin_convolution(_ecm(i), fcn, 0., x_max, efficiency());
    }
  }
  if (isEnergySpreadEnabled()) {
    _integralOperatorMatrix =  _energySpreadMatrix() * _integralOperatorMatrix;
  }
}

TF1* ISRSolverSLAE::_createInterpFunction() const {
  std::function<double(double*, double*)> fcn = [this](double* x, double* par) {
    return this->_interp.eval(this->bcs(), x[0]);
  };
  auto f1 = new TF1("interpFCN", fcn, _energyThreshold(), _ecm(_getN() - 1), 0);
  f1->SetNpx(1.e+4);
  return f1;
}

TF1* ISRSolverSLAE::_createDerivativeInterpFunction(
    std::size_t p, const std::string& name) const {
  std::function<double(double*, double*)> fcn = [p, this](double* x,
                                                          double* par) {
    return this->_interp.derivEval(this->bcs(), x[0]);
  };
  auto f1 =
      new TF1(name.c_str(), fcn, _energyThreshold(), _ecm(_getN() - 1), 0);
  f1->SetNpx(1.e+4);
  return f1;
}

Eigen::VectorXd ISRSolverSLAE::_bcsErr() const {
  return _covMatrixBornCS.diagonal().array().pow(0.5);
}

void ISRSolverSLAE::_evalDotProductOperator() {
  std::size_t i;
  _dotProdOp = Eigen::RowVectorXd(_getN());
  std::function<double(double)> fcn =
      [&i, this](double energy) {
        return this->_interp.basisEval(i, energy);
      };
  for (i = 0; i < _getN(); ++i) {
    // TO-DO : optimize
    double error;
    _dotProdOp(i) = integrate(fcn, getThresholdEnergy(), getMaxEnergy(), error);
  }
}

const Eigen::RowVectorXd& ISRSolverSLAE::_getDotProdOp() const {
  return _dotProdOp;
}

Eigen::MatrixXd ISRSolverSLAE::_energySpreadMatrix() const {
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(_getN(), _getN());
  std::size_t i;
  std::size_t j;
  std::function<double(double)> fcn =
      [&j, this](double energy) {
        double result = this->_interp.basisEval(j, energy);
        return result;
      };
  for (i = 0; i < _getN(); ++i) {
    for (j = 0; j < _getN(); ++j) {
      double sigma2 = std::pow(this->_ecmErr(i), 2);
      result(i, j) = gaussian_conv(this->_ecm(i), sigma2, fcn);
    }
  }
  return result;
}
