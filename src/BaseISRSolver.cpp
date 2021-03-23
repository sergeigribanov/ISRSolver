#include "BaseISRSolver.hpp"

#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TH2F.h>

#include <algorithm>
#include <vector>

BaseISRSolver::BaseISRSolver(const std::string& inputPath,
                             const InputOptions& inputOpts)
  : _energySpread(false),
  _energyT(inputOpts.thresholdEnergy),
    _efficiency([](double, double) {return 1.;}),
    _tefficiency(nullptr) {
  auto fl = TFile::Open(inputPath.c_str(), "read");
  auto graph = dynamic_cast<TGraphErrors*>(
      fl->Get(inputOpts.visibleCSGraphName.c_str()));
  _n = graph->GetN();
  double energyKoeff;
  if (inputOpts.energyUnitMeVs)
    energyKoeff = 1.e-3;
  else
    energyKoeff = 1;
  std::vector<CSData> visibleCS;
  visibleCS.reserve(_n);
  for (std::size_t i = 0; i < _n; ++i)
    visibleCS.push_back({.cmEnergy = graph->GetX()[i] * energyKoeff,
                          .cs = graph->GetY()[i],
                          .cmEnergyError = graph->GetEX()[i] * energyKoeff,
                          .csError = graph->GetEY()[i]});
  if (inputOpts.efficiencyName.length() > 0) {
    _tefficiency = dynamic_cast<TEfficiency*>(fl->Get(inputOpts.efficiencyName.c_str())->Clone("tefficiency"));
  }
  fl->Close();
  delete fl;
  std::sort(
      visibleCS.begin(), visibleCS.end(),
      [](const CSData& x, const CSData& y) { return x.cmEnergy < y.cmEnergy; });
  _visibleCSData = {.cmEnergy = Eigen::VectorXd(_n),
                     .cmEnergyError = Eigen::VectorXd(_n),
                     .cs = Eigen::VectorXd(_n),
                     .csError = Eigen::VectorXd(_n)};
  std::transform(visibleCS.begin(), visibleCS.end(),
                 _visibleCSData.cmEnergy.data(),
                 [](const CSData& x) { return x.cmEnergy; });
  std::transform(visibleCS.begin(), visibleCS.end(),
                 _visibleCSData.cmEnergyError.data(),
                 [](const CSData& x) { return x.cmEnergyError; });
  std::transform(visibleCS.begin(), visibleCS.end(),
                 _visibleCSData.cs.data(),
                 [](const CSData& x) { return x.cs; });
  std::transform(visibleCS.begin(), visibleCS.end(),
                 _visibleCSData.csError.data(),
                 [](const CSData& x) { return x.csError; });
  if (_tefficiency)
    _setupEfficiency();
}

BaseISRSolver::BaseISRSolver(const BaseISRSolver& solver) :
  _energySpread(solver._energySpread),
  _energyT(solver._energyT), _n(solver._n),
  _visibleCSData(solver._visibleCSData),
  _efficiency(solver._efficiency),
  _tefficiency(solver._tefficiency),
  _bornCS(solver._bornCS) {}

BaseISRSolver::~BaseISRSolver() {
  if (_tefficiency) {
    delete _tefficiency;
  }
}

void BaseISRSolver::_setupEfficiency() {
  if (_tefficiency->GetDimension() == 1) {
    _efficiency = [this](double, double energy) {
      int bin = this->_tefficiency->FindFixBin(energy);
      return this->_tefficiency->GetEfficiency(bin);
    };
  }
  if (_tefficiency->GetDimension() == 2) {
    _efficiency = [this](double x, double energy) {
      int bin = this->_tefficiency->FindFixBin(x, energy);
      return this->_tefficiency->GetEfficiency(bin);
    };
  }
  if (_tefficiency->GetDimension() > 2) {
    EfficiencyDimensionException ex;
    throw ex;
  }
}

std::size_t BaseISRSolver::_getN() const { return _n; }

double BaseISRSolver::_energyThreshold() const { return _energyT; }

double BaseISRSolver::_sThreshold() const { return _energyT * _energyT; }

double BaseISRSolver::_s(std::size_t i) const {
  return _visibleCSData.cmEnergy(i) * _visibleCSData.cmEnergy(i);
}

double BaseISRSolver::_ecm(std::size_t i) const {
  return _visibleCSData.cmEnergy(i);
}

const Eigen::VectorXd& BaseISRSolver::ecm() const {
  return _visibleCSData.cmEnergy;
}

Eigen::VectorXd& BaseISRSolver::_ecm() { return _visibleCSData.cmEnergy; }

const Eigen::VectorXd& BaseISRSolver::ecmErr() const {
  return _visibleCSData.cmEnergyError;
}

Eigen::VectorXd& BaseISRSolver::_ecmErr() {
  return _visibleCSData.cmEnergyError;
}

double BaseISRSolver::_ecmErr(std::size_t i) const {
  return  _visibleCSData.cmEnergyError(i);
}

const Eigen::VectorXd& BaseISRSolver::_vcs() const {
  return _visibleCSData.cs;
}

Eigen::VectorXd& BaseISRSolver::_vcs() { return _visibleCSData.cs; }

const Eigen::VectorXd& BaseISRSolver::_vcsErr() const {
  return _visibleCSData.csError;
}

Eigen::VectorXd& BaseISRSolver::_vcsErr() { return _visibleCSData.csError; }

Eigen::MatrixXd BaseISRSolver::_vcsInvCovMatrix() const {
  return _visibleCSData.csError.array().pow(-2.).matrix().asDiagonal();
}

const Eigen::VectorXd& BaseISRSolver::bcs() const { return _bornCS; }

Eigen::VectorXd& BaseISRSolver::_bcs() { return _bornCS; }

const std::function<double(double, double)>& BaseISRSolver::efficiency() const {
  return _efficiency;
}

double BaseISRSolver::getThresholdEnergy() const {
  return _energyT;
}

double BaseISRSolver::getMinEnergy() const {
  return _visibleCSData.cmEnergy(0);
}

double BaseISRSolver::getMaxEnergy() const {
  return _visibleCSData.cmEnergy(_n - 1);
}

bool BaseISRSolver::isEnergySpreadEnabled() const {
  return _energySpread;
}

void BaseISRSolver::enableEnergySpread() {
  _energySpread = true;
}

void BaseISRSolver::disableEnergySpread() {
  _energySpread = false;
}

void BaseISRSolver::resetVisibleCS(const Eigen::VectorXd& vcs) {
  _visibleCSData.cs = vcs;
}

void BaseISRSolver::resetVisibleCSErrors(const Eigen::VectorXd& vcsErr) {
  _visibleCSData.csError = vcsErr;
}

void BaseISRSolver::resetECMErrors(const Eigen::VectorXd& ecmErr) {
  _visibleCSData.cmEnergyError = ecmErr;
}

const Eigen::VectorXd& BaseISRSolver::vcs() const {
  return _visibleCSData.cs;
}

const Eigen::VectorXd& BaseISRSolver::vcsErr() const {
  return _visibleCSData.csError;
}
