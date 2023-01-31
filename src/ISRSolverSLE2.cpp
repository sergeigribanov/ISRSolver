#include "ISRSolverSLE2.hpp"

#include <TFile.h>
#include <TGraphErrors.h>
#include <TMatrixD.h>

#include <algorithm>
#include <boost/lexical_cast.hpp>
#include <cmath>
#include <fstream>
#include <set>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/SVD>

#include "Integration.hpp"
#include "KuraevFadin.hpp"

double* extractIntOpMatrix(ISRSolverSLE2* solver) {
  return solver->_integralOperatorMatrix.data();
}

double* extractBCSCovMatrix(ISRSolverSLE2* solver) {
  return solver->_covMatrixBornCS.data();
}

ISRSolverSLE2::ISRSolverSLE2(std::size_t numberOfPoints,
                             double* energy,
                             double* visibleCS,
                             double* energyErr,
                             double* visibleCSErr,
                             double thresholdEnergy,
                             const std::function<double(double, std::size_t)>&
                             efficiency) :
    BaseISRSolver2(numberOfPoints,
                  energy,
                  visibleCS,
                  energyErr,
                  visibleCSErr,
                  thresholdEnergy,
                  efficiency),
    _interp(Interpolator2(ecm(), getThresholdEnergy())),
    _isEqMatrixPrepared(false) {}

ISRSolverSLE2::ISRSolverSLE2(TGraphErrors* vcsGraph,
                             double thresholdEnergy) :
    BaseISRSolver2(vcsGraph, thresholdEnergy),
    _interp(Interpolator2(ecm(), getThresholdEnergy())),
    _isEqMatrixPrepared(false) {}

ISRSolverSLE2::ISRSolverSLE2(TGraphErrors* vcsGraph,
                             const std::vector<TH1D*>& eff,
                             double thresholdEnergy) :
    BaseISRSolver2(vcsGraph, eff, thresholdEnergy),
    _interp(Interpolator2(ecm(), getThresholdEnergy())),
    _isEqMatrixPrepared(false) {}

ISRSolverSLE2::ISRSolverSLE2(const std::string& inputPath,
                             const InputOptions& inputOpts) :
    BaseISRSolver2(inputPath, inputOpts),
    _interp(Interpolator2(ecm(), getThresholdEnergy())),
    _isEqMatrixPrepared(false) {}

ISRSolverSLE2::ISRSolverSLE2(const ISRSolverSLE2& solver) :
  BaseISRSolver2::BaseISRSolver2(solver),
  _interp(solver._interp),
  _isEqMatrixPrepared(solver._isEqMatrixPrepared),
  _integralOperatorMatrix(solver._integralOperatorMatrix),
  _covMatrixBornCS(solver._covMatrixBornCS),
  _dotProdOp(solver._dotProdOp) {}

ISRSolverSLE2::~ISRSolverSLE2() {}

const Eigen::MatrixXd& ISRSolverSLE2::getIntegralOperatorMatrix() const {
  return _integralOperatorMatrix;
}

const Eigen::MatrixXd& ISRSolverSLE2::getBornCSCovMatrix() const {
  return _covMatrixBornCS;
}

Eigen::MatrixXd& ISRSolverSLE2::_getIntegralOperatorMatrix() {
  return _integralOperatorMatrix;
}

Eigen::MatrixXd& ISRSolverSLE2::_getBornCSCovMatrix() {
  return _covMatrixBornCS;
}

void ISRSolverSLE2::solve() {
  if (!_isEqMatrixPrepared) {
    evalEqMatrix();
    _isEqMatrixPrepared = true;
  }
  _bcs() =
      _integralOperatorMatrix.completeOrthogonalDecomposition().solve(_vcs());
  _covMatrixBornCS = (_integralOperatorMatrix.transpose() *
                      _vcsInvCovMatrix() * _integralOperatorMatrix).inverse();
}

void ISRSolverSLE2::save(const std::string& outputPath,
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

void ISRSolverSLE2::setRangeInterpSettings(
    const std::vector<std::tuple<bool, int, int>>& interpRangeSettings) {
  _interp = Interpolator2(interpRangeSettings, ecm(), getThresholdEnergy());
}

void ISRSolverSLE2::setRangeInterpSettingsJSON(const json& interpRangeSettings) {
  _interp = Interpolator2(interpRangeSettings, ecm(), getThresholdEnergy());
}

void ISRSolverSLE2::setRangeInterpSettings(const std::string& pathToJSON) {
  _interp = Interpolator2(pathToJSON, ecm(), getThresholdEnergy());
}

void ISRSolverSLE2::evalEqMatrix() {
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

TF1* ISRSolverSLE2::_createInterpFunction() const {
  std::function<double(double*, double*)> fcn = [this](double* x, double* par) {
    double result = this->_interp.eval(this->bcs(), x[0]);
    return result;
  };
  auto f1 = new TF1("interpFCN", fcn, _energyThreshold(), _ecm(_getN() - 1), 0);
  f1->SetNpx(1.e+4);
  return f1;
}

Eigen::VectorXd ISRSolverSLE2::_bcsErr() const {
  return _covMatrixBornCS.diagonal().array().pow(0.5);
}

void ISRSolverSLE2::_evalDotProductOperator() {
  std::size_t i;
  _dotProdOp = Eigen::RowVectorXd(_getN());
  for (i = 0; i < _getN(); ++i) {
    _dotProdOp(i) = _interp.evalIntegralBasis(i);
  }
}

const Eigen::RowVectorXd& ISRSolverSLE2::_getDotProdOp() const {
  return _dotProdOp;
}

Eigen::MatrixXd ISRSolverSLE2::_energySpreadMatrix() const {
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

double ISRSolverSLE2::interpEval(const Eigen::VectorXd& y,
                                 double energy) const {
  return _interp.eval(y, energy);
}

double ISRSolverSLE2::evalConditionNumber() const {
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(getIntegralOperatorMatrix());
  return svd.singularValues()(0) /
      svd.singularValues()(svd.singularValues().size()-1);
}

void ISRSolverSLE2::printConditionNumber() const {
  std::cout << "condition number = " << evalConditionNumber() << std::endl;
}

double ISRSolverSLE2::sConvolution(
    const std::function<double(double)>& fcn) const {
  std::function<double(double)> ifcn =
      [fcn, this](double s) {
        const double en = std::sqrt(s);
        const double result = this->_interp.eval(this->bcs(), en) * fcn(s);
        return result;
      };
  double error;
  const double s_min = getThresholdEnergy() * getThresholdEnergy();
  const double s_max = getMaxEnergy() * getMaxEnergy();
  return integrate(ifcn, s_min, s_max, error);
}

double ISRSolverSLE2::sConvolution(
    const std::function<double(double)>& fcn,
    double s_min, double s_max) const {
  std::function<double(double)> ifcn =
      [fcn, this](double s) {
        const double en = std::sqrt(s);
        const double result = this->_interp.eval(this->bcs(), en) * fcn(s);
        return result;
      };
  double error;
  const double s1_min = std::max(getThresholdEnergy() * getThresholdEnergy(), s_min);
  const double s1_max = std::min(getMaxEnergy() * getMaxEnergy(), s_max);
  return integrate(ifcn, s1_min, s1_max, error);
}


Eigen::RowVectorXd ISRSolverSLE2::sConvolutionOperator(
    const std::function<double(double)>& fcn) const {
  Eigen::RowVectorXd result = Eigen::RowVectorXd::Zero(_getN());
  for (std::size_t j = 0; j < _getN(); ++j) {
    result(j) = _interp.evalBasisSConvolution(j, fcn);
  }
  return result;
}

Eigen::RowVectorXd ISRSolverSLE2::sConvolutionOperator(
    const std::function<double(double)>& fcn,
    double s_min, double s_max) const {
  Eigen::RowVectorXd result = Eigen::RowVectorXd::Zero(_getN());
  for (std::size_t j = 0; j < _getN(); ++j) {
    result(j) = _interp.evalBasisSConvolution(j, fcn, s_min, s_max);
  }
  return result;
}
