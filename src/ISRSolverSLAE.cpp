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

#include "integration.hpp"
#include "kuraev_fadin.hpp"

using json = nlohmann::json;

ISRSolverSLAE::ISRSolverSLAE(const std::string& inputPath,
                             const InputOptions& inputOpts) :
  BaseISRSolver(inputPath, inputOpts) {
  _setDefaultInterpSettings();
}

ISRSolverSLAE::ISRSolverSLAE(const ISRSolverSLAE& solver) :
  BaseISRSolver::BaseISRSolver(solver),
  _interpSettings(solver._interpSettings),
  _integralOperatorMatrix(solver._integralOperatorMatrix),
  _invBornCSErrorMatrix(solver._invBornCSErrorMatrix) {
}

ISRSolverSLAE::~ISRSolverSLAE() {
}

const Eigen::MatrixXd& ISRSolverSLAE::getIntegralOperatorMatrix() const {
  return _integralOperatorMatrix;
}

const Eigen::MatrixXd& ISRSolverSLAE::getInverseBornCSErrorMatrix() const {
  return _invBornCSErrorMatrix;
}

Eigen::MatrixXd& ISRSolverSLAE::_getIntegralOperatorMatrix() {
  return _integralOperatorMatrix;
}

Eigen::MatrixXd& ISRSolverSLAE::_getInverseBornCSErrorMatrix() {
  return _invBornCSErrorMatrix;
}

void ISRSolverSLAE::solve() {
  _evalEqMatrix();
  _bcs() =
      _integralOperatorMatrix.completeOrthogonalDecomposition().solve(_vcs());
  _invBornCSErrorMatrix = _integralOperatorMatrix.transpose() *
                          _vcsInvErrMatrix() * _integralOperatorMatrix;
}

void ISRSolverSLAE::save(const std::string& outputPath,
                         const OutputOptions& outputOpts) {
  TGraphErrors vcs(_getN(), _ecm().data(), _vcs().data(), _ecmErr().data(),
                   _vcsErr().data());
  TGraphErrors bcs(_getN(), _ecm().data(), _bcs().data(), _ecmErr().data(),
                   _bcsErr().data());
  TMatrixD intergalOperatorMatrix(_getN(), _getN());
  Eigen::MatrixXd tmpIntOpM = _integralOperatorMatrix.transpose();
  intergalOperatorMatrix.SetMatrixArray(tmpIntOpM.data());
  TMatrixD bornCSInverseErrorMatrix(_getN(), _getN());
  Eigen::MatrixXd tmpInvErrM = _invBornCSErrorMatrix.transpose();
  bornCSInverseErrorMatrix.SetMatrixArray(tmpInvErrM.data());
  auto f0 = _createInterpFunction();
  auto f1 = _createDerivativeInterpFunction(1, "interp1DivFCN");
  auto fl = TFile::Open(outputPath.c_str(), "recreate");
  fl->cd();
  vcs.Write(outputOpts.visibleCSGraphName.c_str());
  bcs.Write(outputOpts.bornCSGraphName.c_str());
  intergalOperatorMatrix.Write("intergalOperatorMatrix");
  bornCSInverseErrorMatrix.Write("bornCSInverseErrorMatrix");
  f0->Write();
  f1->Write();
  fl->Close();
  delete f0;
  delete f1;
  delete fl;
}

void ISRSolverSLAE::setInterpSettings(
    const std::vector<InterpPtSettings>& interpSettings) {
  if (interpSettings.size() != _getN()) {
    InterpSettingsSizeException ex;
    throw ex;
  }
  _interpSettings = interpSettings;
}

void ISRSolverSLAE::setInterpSettings(const std::string& pathToJSON) {
  std::ifstream fl(pathToJSON);
  json s;
  fl >> s;
  std::vector<InterpPtSettings> interpSettings(_getN());
  std::set<std::size_t> keys;
  std::size_t key;
  for (const auto& el : s.items()) {
    key = boost::lexical_cast<std::size_t>(el.key());
    keys.insert(key);
    interpSettings[key] = {.numCoeffs = el.value()[0],
                           .nPointsLeft = el.value()[1]};
  }
  if (keys.size() != _getN()) {
    InterpSettingsSizeException ex;
    throw ex;
  }
  _interpSettings = interpSettings;
}

double ISRSolverSLAE::_getXmin(int i, int j) const { return 1 - _s(j) / _s(i); }

double ISRSolverSLAE::_getXmax(int i, int j) const {
  if (j == 0) return 1 - _sThreshold() / _s(i);
  return 1 - _s(j - 1) / _s(i);
}

double ISRSolverSLAE::_getNumCoeffs(int j) const {
  return _interpSettings[j].numCoeffs;
}

double ISRSolverSLAE::_getNPointsLeft(int j) const {
  return _interpSettings[j].nPointsLeft;
}

Eigen::MatrixXd ISRSolverSLAE::_permutation(int j) const {
  const std::size_t nc = _getNumCoeffs(j);
  const std::size_t npl = _getNPointsLeft(j);
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(nc, _getN());
  int p = j - npl;
  for (std::size_t l = 0; l < nc; ++l) {
    if (p >= 0) {
      result(l, p) = 1;
    }
    p++;
  }
  return result;
}

Eigen::MatrixXd ISRSolverSLAE::_interpInvMatrix(int j) const {
  const std::size_t nc = _getNumCoeffs(j);
  const std::size_t npl = _getNPointsLeft(j);
  Eigen::MatrixXd interpMatrix = Eigen::MatrixXd::Zero(nc, nc);
  int p = j - npl;
  for (std::size_t l = 0; l < nc; ++l) {
    if (p < 0) {
      for (std::size_t k = 0; k < nc; ++k)
        interpMatrix(l, k) = std::pow(_sThreshold(), k);
    } else {
      for (std::size_t k = 0; k < nc; ++k)
        interpMatrix(l, k) = std::pow(_s(p), k);
    }
    p++;
  }
  return interpMatrix.inverse();
}

Eigen::MatrixXd ISRSolverSLAE::_polConvKuraevFadinMatrix(int j) const {
  const std::size_t nc = _getNumCoeffs(j);
  Eigen::MatrixXd result(_getN(), nc);
  std::size_t k;
  std::function<double(double)> fcn = [&k, this](double t) {return std::pow(t, k); };
  for (std::size_t i = j; i < _getN(); ++i) {
    for (k = 0; k < nc; ++k) {
      result(i, k) =
	kuraev_fadin_convolution(_s(i), fcn, _getXmin(i, j), _getXmax(i, j), efficiency());
    }
  }
  return result;
}

Eigen::MatrixXd ISRSolverSLAE::_evalA(int j) const {
  return _polConvKuraevFadinMatrix(j) * _interpInvMatrix(j) * _permutation(j);
}

void ISRSolverSLAE::_evalEqMatrix() {
  std::vector<Eigen::MatrixXd> ai;
  ai.reserve(_getN());
  Eigen::MatrixXd tmpV = Eigen::MatrixXd::Zero(_getN(), _getN());
  for (std::size_t i = 0; i < _getN(); ++i) {
    tmpV += _evalA(i);
    ai.push_back(tmpV);
  }
  _integralOperatorMatrix = Eigen::MatrixXd::Zero(_getN(), _getN());
  for (std::size_t i = 0; i < _getN(); ++i) {
    for (std::size_t p = 0; p < _getN(); ++p) {
      _integralOperatorMatrix(i, p) += ai[i](i, p);
    }
  }
}

TF1* ISRSolverSLAE::_createInterpFunction() const {
  std::function<double(double*, double*)> fcn = [this](double* x, double* par) {
    double en = x[0];
    return this->_interpProjector(en) * this->bcs();
  };
  auto f1 = new TF1("interpFCN", fcn, _energyThreshold(), _ecm(_getN() - 1), 0);
  f1->SetNpx(1.e+4);
  return f1;
}

TF1* ISRSolverSLAE::_createDerivativeInterpFunction(
    std::size_t p, const std::string& name) const {
  std::function<double(double*, double*)> fcn = [p, this](double* x,
                                                          double* par) {
    double en = x[0];
    return this->_interpDerivativeProjector(en, p) * this->bcs();
  };
  auto f1 =
      new TF1(name.c_str(), fcn, _energyThreshold(), _ecm(_getN() - 1), 0);
  f1->SetNpx(1.e+4);
  return f1;
}

Eigen::RowVectorXd ISRSolverSLAE::_interpProjector(double en) const {
  Eigen::RowVectorXd result = Eigen::RowVectorXd::Zero(_getN());
  if (en <= _energyThreshold()) return result;
  std::size_t i = 0;
  double enp = _energyThreshold();
  while (i < _getN()) {
    if (en > enp && en <= _ecm(i)) break;
    enp = _ecm(i);
    i++;
  }
  if (i == _getN()) i--;
  Eigen::MatrixXd coeffs = _interpInvMatrix(i) * _permutation(i);
  const std::size_t nc = _getNumCoeffs(i);
  for (std::size_t k = 0; k < nc; ++k)
    result += coeffs.row(k) * std::pow(en * en, k);
  return result;
}

Eigen::RowVectorXd ISRSolverSLAE::_interpDerivativeProjector(
    double en, std::size_t p) const {
  Eigen::RowVectorXd result = Eigen::RowVectorXd::Zero(_getN());
  if (en <= _energyThreshold()) return result;
  std::size_t i = 0;
  double enp = _energyThreshold();
  while (i < _getN()) {
    if (en > enp && en <= _ecm(i)) break;
    enp = _ecm(i);
    i++;
  }
  if (i == _getN()) i--;
  Eigen::MatrixXd coeffs = _interpInvMatrix(i) * _permutation(i);
  const std::size_t nc = _getNumCoeffs(i);
  for (std::size_t k = p; k < nc; ++k) {
    double v = 1;
    for (std::size_t o = 0; o < p; ++o) v *= k - o;
    result += coeffs.row(k) * v * std::pow(en * en, k - p);
  }
  return result;
}

Eigen::VectorXd ISRSolverSLAE::_bcsErr() const {
  return _invBornCSErrorMatrix.diagonal().array().pow(-0.5);
}

void ISRSolverSLAE::_setDefaultInterpSettings() {
  _interpSettings.resize(_getN(), {.numCoeffs = 2, .nPointsLeft = 1});
}
