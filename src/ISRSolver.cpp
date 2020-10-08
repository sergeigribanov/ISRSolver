#include "ISRSolver.hpp"

#include <TFile.h>
#include <TGraphErrors.h>
#include <TMatrixD.h>

#include <algorithm>
#include <boost/lexical_cast.hpp>
#include <cmath>
#include <fstream>
#include <functional>
#include <nlohmann/json.hpp>
#include <set>
#include <nlopt.hpp>

#include "kuraev_fadin.hpp"
#include "integration.hpp"

using json = nlohmann::json;

double objectiveFCN(unsigned n, const double* z, double* grad, void* solver) {
  auto sp = reinterpret_cast<ISRSolver*>(solver);
  Eigen::Map<const Eigen::VectorXd> vm(z, n);
  if (grad) {
    Eigen::VectorXd vgrad = sp->evalRegFuncGradNorm2(vm);
    for (std::size_t i = 0; i < n; ++i) {
      grad[i] = vgrad(i);
    }
  }
  return sp->evalRegFuncNorm2(vm);
}

ISRSolver::ISRSolver(const std::string& inputPath,
                     const InputOptions& inputOpts)
  : _alpha(1.e-4), _beta(_alpha), _inputOpts(inputOpts),
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
    return this->interpProjector(en) * this->_bornCS;
  };
  auto f1 = new TF1("interpFCN", fcn, _inputOpts.thresholdEnergy,
                    _measuredCSData.cmEnergy(getN() - 1), 0);
  f1->SetNpx(1.e+4);
  return f1;
}

TF1* ISRSolver::createDerivativeInterpFunction(std::size_t p, const std::string& name) const {
  std::function<double(double*, double*)> fcn = [p, this](double* x, double* par) {
    double en = x[0];
    return this->interpDerivativeProjector(en, p) * this->_bornCS;
  };
  auto f1 = new TF1(name.c_str(), fcn, _inputOpts.thresholdEnergy,
                    _measuredCSData.cmEnergy(getN() - 1), 0);
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

void ISRSolver::setInterpSettings(const std::string& pathToJSON) {
  std::ifstream fl(pathToJSON);
  json s;
  fl >> s;
  std::vector<InterpPtSettings> interpSettings(getN());
  std::set<std::size_t> keys;
  std::size_t key;
  for (const auto& el : s.items()) {
    key = boost::lexical_cast<std::size_t>(el.key());
    keys.insert(key);
    interpSettings[key] = {.numCoeffs = el.value()[0],
                           .nPointsLeft = el.value()[1]};
  }
  if (keys.size() != getN()) {
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
  auto f0 = createInterpFunction();
  auto f1 = createDerivativeInterpFunction(1, "interp1DivFCN");
  Eigen::VectorXd vIntegral = evalIntegralMatrix() * _bornCS;
  TGraphErrors gIntegral(getN(), _measuredCSData.cmEnergy.data(), vIntegral.data(), 0, 0);
  Eigen::VectorXd vDerivative = interpPointWiseDerivativeProjector() * _bornCS;
  TGraphErrors gDerivative(getN(), _measuredCSData.cmEnergy.data(), vDerivative.data(), 0, 0);
  
  auto fl = TFile::Open(outputPath.c_str(), "recreate");
  fl->cd();
  vcs.Write(outputOpts.measuredCSGraphName.c_str());
  bcs.Write(outputOpts.bornCSGraphName.c_str());
  intergalOperatorMatrix.Write("intergalOperatorMatrix");
  bornCSInverseErrorMatrix.Write("bornCSInverseErrorMatrix");
  gIntegral.Write("bcs_integral");
  gDerivative.Write("bcs_deriv");
  f0->Write();
  f1->Write();
  fl->Close();
  delete f0;
  delete f1;
  delete fl;
}

Eigen::RowVectorXd ISRSolver::interpProjector(double en) const {
  Eigen::RowVectorXd result = Eigen::RowVectorXd::Zero(getN());
  if (en <= _inputOpts.thresholdEnergy) return result;
  std::size_t i = 0;
  double enp = _inputOpts.thresholdEnergy;
  while (i < getN()) {
    if (en > enp && en <= _measuredCSData.cmEnergy(i)) break;
    enp = _measuredCSData.cmEnergy(i);
    i++;
  }
  if (i == getN()) i--;
  Eigen::MatrixXd coeffs = interpInvMatrix(i) * permutation(i);
  const std::size_t nc = getNumCoeffs(i);
  for (std::size_t k = 0; k < nc; ++k)
    result += coeffs.row(k) * std::pow(en * en, k);
  return result;
}

Eigen::RowVectorXd ISRSolver::interpDerivativeProjector(double en, std::size_t p) const {
  Eigen::RowVectorXd result = Eigen::RowVectorXd::Zero(getN());
  if (en <= _inputOpts.thresholdEnergy) return result;
  std::size_t i = 0;
  double enp = _inputOpts.thresholdEnergy;
  while (i < getN()) {
    if (en > enp && en <= _measuredCSData.cmEnergy(i)) break;
    enp = _measuredCSData.cmEnergy(i);
    i++;
  }
  if (i == getN()) i--;
  Eigen::MatrixXd coeffs = interpInvMatrix(i) * permutation(i);
  const std::size_t nc = getNumCoeffs(i);
  for (std::size_t k = p; k < nc; ++k) {
    double v = 1;
    for (std::size_t o = 0; o < p; ++o) v *= k - o;
    result += coeffs.row(k) * v * std::pow(en * en, k - p);
  }
  return result;
}

Eigen::MatrixXd ISRSolver::interpPointWiseDerivativeProjector() const {
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(getN(), getN());
  for (std::size_t j = 0; j < getN(); ++j) {
      Eigen::MatrixXd coeffs = interpInvMatrix(j) * permutation(j);
      const std::size_t nc = getNumCoeffs(j);
      for (std::size_t k = 1; k < nc; ++k) {
	result.row(j) += coeffs.row(k) * k * std::pow(_measuredCSData.s(j), k - 1);
      }
  }
  return result;
}

Eigen::RowVectorXd ISRSolver::polIntegralOp(int j) const {
  const std::size_t nc = getNumCoeffs(j);
  Eigen::VectorXd result(nc);
  std::size_t k;
  double error;
  double sMin;
  if (j == 0) {
    sMin = _sT;
  } else {
    sMin = _measuredCSData.s(j - 1);
  }
  double sMax = _measuredCSData.s(j);
  std::function<double(double)> fcn = [&k](double t) { return std::pow(t, k); };
  for (k = 0; k < nc; ++k) {
    result(k) = integrate(fcn, sMin, sMax, error);
  }
  return result.transpose();
}

Eigen::RowVectorXd ISRSolver::evalPolA(int j) const {
  return polIntegralOp(j) * interpInvMatrix(j) * permutation(j);
}

Eigen::MatrixXd ISRSolver::evalIntegralMatrix() const {
  std::vector<Eigen::RowVectorXd> ai;
  ai.reserve(getN());
  Eigen::MatrixXd tmpV = Eigen::RowVectorXd::Zero(getN());
  for (std::size_t i = 0; i < getN(); ++i) {
    tmpV += evalPolA(i);
    ai.push_back(tmpV);
  }
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(getN(), getN());
  for (std::size_t i = 0; i < getN(); ++i) {
    for (std::size_t p = 0; p < getN(); ++p) {
      result(i, p) += ai[i](p);
    }
  }
  return result;
}

Eigen::RowVectorXd ISRSolver::evalScalarProductOperator() const {
  Eigen::MatrixXd result = Eigen::RowVectorXd::Zero(getN());
  for (std::size_t i = 0; i < getN(); ++i) {
    result += evalPolA(i);
  }
  return result;
}

double ISRSolver::evalRegFuncNorm2(const Eigen::VectorXd& z) const {
  // TO DO: add derivative, insert integral norm
  Eigen::RowVectorXd op = evalScalarProductOperator();
  Eigen::VectorXd dv = _integralOperatorMatrix * z - _measuredCSData.cs;
  Eigen::VectorXd dz = interpPointWiseDerivativeProjector() * z;
  return op * (dv.array() * dv.array() + _alpha * z.array() * z.array() +
	       _beta * dz.array() * dz.array()).matrix();
}

Eigen::VectorXd ISRSolver::evalRegFuncGradNorm2(const Eigen::VectorXd& z) const {
  Eigen::VectorXd dv = _integralOperatorMatrix * z - _measuredCSData.cs;
  Eigen::RowVectorXd op = evalScalarProductOperator();
  Eigen::VectorXd result = Eigen::VectorXd::Zero(getN());
  Eigen::MatrixXd dzOp = interpPointWiseDerivativeProjector();
  Eigen::VectorXd dz = dzOp * z;
  for (std::size_t l = 0; l < getN(); ++l) {
    for (std::size_t i = 0; i < getN(); ++i) {
      result(l) += 2 * op(i) * _integralOperatorMatrix(i, l) * dv(i) +
	2 * _beta * op(i) * dzOp(i, l) * dz(i);
      if (l == i) {
	result(l) += 2 * _alpha * op(l) * z(l);
      }
    }
  }
  return result;
}

void ISRSolver::solveTikhonov() {
  _integralOperatorMatrix = evalEqMatrix();
  std::vector<double> lowerBounds(_n, 0);
  nlopt::opt opt(nlopt::LD_MMA, _n);
  opt.set_lower_bounds(lowerBounds);
  opt.set_min_objective(objectiveFCN, this);
  opt.set_xtol_rel(1.e-6);
  std::vector<double> z(_n, 0);
  double minf;
  try {
    opt.optimize(z, minf);
    _bornCS = Eigen::Map<Eigen::VectorXd>(z.data(), _n);
    std::cout << _bornCS << std::endl;
    TGraphErrors bcs(getN(), _measuredCSData.cmEnergy.data(), _bornCS.data(),
		     _measuredCSData.cmEnergyError.data(), 0);
    auto fl = TFile::Open("bcs_test.root", "recreate");
    fl->cd();
    bcs.Write("bcs");
    fl->Close();
  }
  catch (std::exception &e) {
    std::cout << "nlopt failed: " << e.what() << std::endl;
  }
}
