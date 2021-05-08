#include "ISRSolverSLE.hpp"

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
#include <Eigen/Core>
#include <Eigen/SVD>

#include "Integration.hpp"
#include "KuraevFadin.hpp"

using json = nlohmann::json;

ISRSolverSLE::ISRSolverSLE(std::size_t numberOfPoints,
                           double* energy, double* visibleCS,
                           double* energyErr, double* visibleCSErr,
                           double thresholdEnergy,
                           const std::function<double(double, double)>&
                           efficiency) :
    BaseISRSolver(numberOfPoints,
                  energy, visibleCS,
                  energyErr, visibleCSErr,
                  thresholdEnergy,
                  efficiency),
    _interp(Interpolator(ecm(), getThresholdEnergy())),
    _isEqMatrixPrepared(false) {}

ISRSolverSLE::ISRSolverSLE(TGraphErrors* vcsGraph,
                           double thresholdEnergy) :
    BaseISRSolver(vcsGraph, thresholdEnergy),
    _interp(Interpolator(ecm(), getThresholdEnergy())),
    _isEqMatrixPrepared(false) {}

ISRSolverSLE::ISRSolverSLE(TGraphErrors* vcsGraph,
                           TEfficiency* eff,
                           double thresholdEnergy) :
    BaseISRSolver(vcsGraph, eff, thresholdEnergy),
    _interp(Interpolator(ecm(), getThresholdEnergy())),
    _isEqMatrixPrepared(false) {}

ISRSolverSLE::ISRSolverSLE(const std::string& inputPath,
                             const InputOptions& inputOpts) :
    BaseISRSolver(inputPath, inputOpts),
    _interp(Interpolator(ecm(), getThresholdEnergy())),
    _isEqMatrixPrepared(false) {}

ISRSolverSLE::ISRSolverSLE(const ISRSolverSLE& solver) :
  BaseISRSolver::BaseISRSolver(solver),
  _interp(solver._interp),
  _isEqMatrixPrepared(solver._isEqMatrixPrepared),
  _integralOperatorMatrix(solver._integralOperatorMatrix),
  _covMatrixBornCS(solver._covMatrixBornCS),
  _dotProdOp(solver._dotProdOp) {}

ISRSolverSLE::~ISRSolverSLE() {}

const Eigen::MatrixXd& ISRSolverSLE::getIntegralOperatorMatrix() const {
  return _integralOperatorMatrix;
}

const Eigen::MatrixXd& ISRSolverSLE::getBornCSCovMatrix() const {
  return _covMatrixBornCS;
}

Eigen::MatrixXd& ISRSolverSLE::_getIntegralOperatorMatrix() {
  return _integralOperatorMatrix;
}

Eigen::MatrixXd& ISRSolverSLE::_getBornCSCovMatrix() {
  return _covMatrixBornCS;
}

void ISRSolverSLE::solve() {
  if (!_isEqMatrixPrepared) {
    evalEqMatrix();
    _isEqMatrixPrepared = true;
  }
  _bcs() =
      _integralOperatorMatrix.completeOrthogonalDecomposition().solve(_vcs());
  _covMatrixBornCS = (_integralOperatorMatrix.transpose() *
                      _vcsInvCovMatrix() * _integralOperatorMatrix).inverse();
}

void ISRSolverSLE::save(const std::string& outputPath,
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
  TMatrixD bornCSInvCovMatrix(_getN(), _getN());
  Eigen::MatrixXd tmpInvCovM = _covMatrixBornCS.inverse().transpose();
  bornCSInvCovMatrix.SetMatrixArray(tmpInvCovM.data());
  auto f0 = _createInterpFunction();
  auto fl = TFile::Open(outputPath.c_str(), "recreate");
  fl->cd();
  vcs.Write(outputOpts.visibleCSGraphName.c_str());
  bcs.Write(outputOpts.bornCSGraphName.c_str());
  intergalOperatorMatrix.Write("intergalOperatorMatrix");
  bornCSCovMatrix.Write("covMatrixBornCS");
  bornCSInvCovMatrix.Write("invCovMatrixBornCS");
  f0->Write();
  fl->Close();
  delete f0;
  delete fl;
}

void ISRSolverSLE::setRangeInterpSettings(
    const std::vector<std::tuple<bool, int, int>>& interpRangeSettings) {
  _interp = Interpolator(interpRangeSettings, ecm(), getThresholdEnergy());
}

void ISRSolverSLE::setRangeInterpSettings(const std::string& pathToJSON) {
  _interp = Interpolator(pathToJSON, ecm(), getThresholdEnergy());
}

void ISRSolverSLE::evalEqMatrix() {
  _integralOperatorMatrix = Eigen::MatrixXd::Zero(_getN(), _getN());
  // TO DO: optimize
  for (std::size_t j = 0; j < _getN(); ++j) {
    for(std::size_t i = 0; i < _getN(); ++i) {
      _integralOperatorMatrix(i, j) = _interp.evalKuraevFadinBasisIntegral(i, j, efficiency());
    }
  }
  if (isEnergySpreadEnabled()) {
    _integralOperatorMatrix =  _energySpreadMatrix() * _integralOperatorMatrix;
  }
}

TF1* ISRSolverSLE::_createInterpFunction() const {
  std::function<double(double*, double*)> fcn = [this](double* x, double* par) {
    return this->_interp.eval(this->bcs(), x[0]);
  };
  auto f1 = new TF1("interpFCN", fcn, _energyThreshold(), _ecm(_getN() - 1), 0);
  f1->SetNpx(1.e+4);
  return f1;
}

Eigen::VectorXd ISRSolverSLE::_bcsErr() const {
  return _covMatrixBornCS.diagonal().array().pow(0.5);
}

void ISRSolverSLE::_evalDotProductOperator() {
  std::size_t i;
  _dotProdOp = Eigen::RowVectorXd(_getN());
  for (i = 0; i < _getN(); ++i) {
    _dotProdOp(i) = _interp.evalIntegralBasis(i);
  }
}

const Eigen::RowVectorXd& ISRSolverSLE::_getDotProdOp() const {
  return _dotProdOp;
}

Eigen::MatrixXd ISRSolverSLE::_energySpreadMatrix() const {
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

double ISRSolverSLE::interpEval(const Eigen::VectorXd& y,
                                 double energy) const {
  return _interp.eval(y, energy);
}

double ISRSolverSLE::evalConditionNumber() const {
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(getIntegralOperatorMatrix());
  return svd.singularValues()(0) /
      svd.singularValues()(svd.singularValues().size()-1);
}

void ISRSolverSLE::printConditionNumber() const {
  std::cout << "condition number = " << evalConditionNumber() << std::endl;
}
