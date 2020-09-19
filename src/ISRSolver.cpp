#include "ISRSolver.hpp"

#include <TFile.h>
#include <TGraphErrors.h>
#include <TMatrixD.h>

#include <algorithm>
#include <cmath>
#include <functional>

#include "kuraev_fadin.hpp"

ISRSolver::ISRSolver(const std::string& inputPath,
                     const InputOptions& inputOpts)
    : _inputOpts(inputOpts),
      _sT(inputOpts.thresholdEnergy * inputOpts.thresholdEnergy) {
  auto fl = TFile::Open(inputPath.c_str(), "read");
  auto graph = dynamic_cast<TGraphErrors*>(
      fl->Get(_inputOpts.measuredCSGraphName.c_str()));
  _n = graph->GetN();
  double energyKoeff;
  if (_inputOpts.energyUnitMeVs)
    energyKoeff = 1.e-3;
  else
    energyKoeff = 1;
  std::vector<CSData> measuredCS;
  measuredCS.reserve(_n);
  for (std::size_t i = 0; i < _n; ++i)
    measuredCS.push_back({.cmEnergy = graph->GetX()[i] * energyKoeff,
                          .cs = graph->GetY()[i],
                          .cmEnergyError = graph->GetEX()[i] * energyKoeff,
                          .csError = graph->GetEY()[i]});
  fl->Close();
  delete fl;
  std::sort(
      measuredCS.begin(), measuredCS.end(),
      [](const CSData& x, const CSData& y) { return x.cmEnergy < y.cmEnergy; });
  _measuredCSData = {.cmEnergy = Eigen::VectorXd(_n),
                     .s = Eigen::VectorXd(_n),
                     .cs = Eigen::VectorXd(_n),
                     .cmEnergyError = Eigen::VectorXd(_n),
                     .invCSErrMatrix = Eigen::MatrixXd::Zero(_n, _n)};
  std::transform(measuredCS.begin(), measuredCS.end(),
                 _measuredCSData.cmEnergy.data(),
                 [](const CSData& x) { return x.cmEnergy; });
  std::transform(measuredCS.begin(), measuredCS.end(), _measuredCSData.s.data(),
                 [](const CSData& x) { return x.cmEnergy * x.cmEnergy; });
  std::transform(measuredCS.begin(), measuredCS.end(),
                 _measuredCSData.cs.data(),
                 [](const CSData& x) { return x.cs; });
  std::transform(measuredCS.begin(), measuredCS.end(),
                 _measuredCSData.cmEnergyError.data(),
                 [](const CSData& x) { return x.cmEnergyError; });
  Eigen::VectorXd tmpv(_n);
  std::transform(measuredCS.begin(), measuredCS.end(), tmpv.data(),
                 [](const CSData& x) { return 1. / x.csError / x.csError; });
  _measuredCSData.invCSErrMatrix.diagonal() = tmpv;
  setDefaultInterpSettings();
}

ISRSolver::~ISRSolver() {}

TF1* ISRSolver::createInterpFunction() const {
  std::function<double(double*, double*)> fcn = [this](double* x, double* par) {
    double en = x[0];
    if (en <= _inputOpts.thresholdEnergy) return 0.;
    std::size_t i = 0;
    double enp = _inputOpts.thresholdEnergy;
    while (i < this->getN()) {
      if (en > enp && en <= _measuredCSData.cmEnergy(i)) break;
      enp = _measuredCSData.cmEnergy(i);
      i++;
    }
    if (i == this->getN()) i--;

    Eigen::VectorXd coeffs =
        this->interpInvMatrix(i) * this->permutation(i) * this->_bornCS;
    double result = 0;
    for (int k = 0; k < coeffs.size(); ++k)
      result += coeffs(k) * std::pow(en * en, k);
    return result;
  };
  auto f1 = new TF1("interpFCN", fcn, this->_inputOpts.thresholdEnergy,
                    this->_measuredCSData.cmEnergy(getN() - 1), 0);
  f1->SetNpx(1.e+4);
  return f1;
}

double ISRSolver::getXmin(int i, int j) const {
  return 1 - _measuredCSData.s(j) / _measuredCSData.s(i);
}

double ISRSolver::getXmax(int i, int j) const {
  if (j == 0) return 1 - _sT / _measuredCSData.s(i);
  return 1 - _measuredCSData.s(j - 1) / _measuredCSData.s(i);
}

std::size_t ISRSolver::getN() const { return _n; }

void ISRSolver::setDefaultInterpSettings() {
  _interpSettings.resize(_n, {.numCoeffs = 2, .nPointsLeft = 1});
}

void ISRSolver::setInterpSettings(
    const std::vector<InterpPtSettings>& interpSettings) {
  if (interpSettings.size() != getN()) {
    InterpSettingsSizeException ex;
    throw ex;
  }
  _interpSettings = interpSettings;
}

double ISRSolver::getNumCoeffs(int j) const {
  return _interpSettings[j].numCoeffs;
}

double ISRSolver::getNPointsLeft(int j) const {
  return _interpSettings[j].nPointsLeft;
}

Eigen::MatrixXd ISRSolver::permutation(int j) const {
  const std::size_t nc = getNumCoeffs(j);
  const std::size_t npl = getNPointsLeft(j);
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(nc, getN());
  int p = j - npl;
  for (std::size_t l = 0; l < nc; ++l) {
    if (p >= 0) {
      result(l, p) = 1;
    }
    p++;
  }
  return result;
}

Eigen::MatrixXd ISRSolver::interpInvMatrix(int j) const {
  const std::size_t nc = getNumCoeffs(j);
  const std::size_t npl = getNPointsLeft(j);
  Eigen::MatrixXd interpMatrix = Eigen::MatrixXd::Zero(nc, nc);
  int p = j - npl;
  for (std::size_t l = 0; l < nc; ++l) {
    if (p < 0) {
      for (std::size_t k = 0; k < nc; ++k)
        interpMatrix(l, k) = std::pow(_sT, k);
    } else {
      for (std::size_t k = 0; k < nc; ++k)
        interpMatrix(l, k) = std::pow(_measuredCSData.s(p), k);
    }
    p++;
  }
  return interpMatrix.inverse();
}

Eigen::MatrixXd ISRSolver::polConvKuraevFadinMatrix(int j) const {
  const std::size_t nc = getNumCoeffs(j);
  Eigen::MatrixXd result(getN(), nc);
  std::size_t k;
  std::function<double(double)> fcn = [&k](double t) { return std::pow(t, k); };
  for (std::size_t i = j; i < getN(); ++i) {
    for (k = 0; k < nc; ++k) {
      result(i, k) = kuraev_fadin_convolution(_measuredCSData.s(i), fcn,
                                              getXmin(i, j), getXmax(i, j));
    }
  }
  return result;
}

Eigen::MatrixXd ISRSolver::evalA(int j) const {
  return polConvKuraevFadinMatrix(j) * interpInvMatrix(j) * permutation(j);
}

Eigen::MatrixXd ISRSolver::evalEqMatrix() const {
  std::vector<Eigen::MatrixXd> ai;
  ai.reserve(getN());
  Eigen::MatrixXd tmpV = Eigen::MatrixXd::Zero(getN(), getN());
  for (std::size_t i = 0; i < getN(); ++i) {
    tmpV += evalA(i);
    ai.push_back(tmpV);
  }
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(getN(), getN());
  for (std::size_t i = 0; i < getN(); ++i) {
    for (std::size_t p = 0; p < getN(); ++p) {
      result(i, p) += ai[i](i, p);
    }
  }
  return result;
}

void ISRSolver::solve() {
  _integralOperatorMatrix = evalEqMatrix();
  _bornCS = _integralOperatorMatrix.completeOrthogonalDecomposition().solve(
      _measuredCSData.cs);
  _invBornCSErrorMatrix = _integralOperatorMatrix.transpose() *
                          _measuredCSData.invCSErrMatrix *
                          _integralOperatorMatrix;
}

void ISRSolver::save(const std::string& outputPath,
                     const OutputOptions& outputOpts) {
  Eigen::VectorXd bornCSError =
      _invBornCSErrorMatrix.diagonal().array().pow(-0.5);
  Eigen::VectorXd measuredCSError =
      _measuredCSData.invCSErrMatrix.diagonal().array().pow(-0.5);
  TGraphErrors vcs(
      getN(), _measuredCSData.cmEnergy.data(), _measuredCSData.cs.data(),
      _measuredCSData.cmEnergyError.data(), measuredCSError.data());
  TGraphErrors bcs(getN(), _measuredCSData.cmEnergy.data(), _bornCS.data(),
                   _measuredCSData.cmEnergyError.data(), bornCSError.data());
  TMatrixD intergalOperatorMatrix(getN(), getN());
  Eigen::MatrixXd tmpIntOpM = _integralOperatorMatrix.transpose();
  intergalOperatorMatrix.SetMatrixArray(tmpIntOpM.data());
  TMatrixD bornCSInverseErrorMatrix(getN(), getN());
  Eigen::MatrixXd tmpInvErrM = _invBornCSErrorMatrix.transpose();
  bornCSInverseErrorMatrix.SetMatrixArray(tmpInvErrM.data());
  auto f1 = createInterpFunction();
  auto fl = TFile::Open(outputPath.c_str(), "recreate");
  fl->cd();
  vcs.Write(outputOpts.measuredCSGraphName.c_str());
  bcs.Write(outputOpts.bornCSGraphName.c_str());
  intergalOperatorMatrix.Write("intergalOperatorMatrix");
  bornCSInverseErrorMatrix.Write("bornCSInverseErrorMatrix");
  f1->Write();
  fl->Close();
  delete f1;
  delete fl;
}
