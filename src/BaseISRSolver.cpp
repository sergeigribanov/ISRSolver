#include "BaseISRSolver.hpp"

#include <TFile.h>
#include <TGraphErrors.h>

#include <algorithm>
#include <vector>

BaseISRSolver::BaseISRSolver(const std::string& inputPath,
                             const InputOptions& inputOpts)
    : _energyT(inputOpts.thresholdEnergy) {
  auto fl = TFile::Open(inputPath.c_str(), "read");
  auto graph = dynamic_cast<TGraphErrors*>(
      fl->Get(inputOpts.measuredCSGraphName.c_str()));
  _n = graph->GetN();
  double energyKoeff;
  if (inputOpts.energyUnitMeVs)
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
                     .cmEnergyError = Eigen::VectorXd(_n),
                     .cs = Eigen::VectorXd(_n),
                     .csError = Eigen::VectorXd(_n)};
  std::transform(measuredCS.begin(), measuredCS.end(),
                 _measuredCSData.cmEnergy.data(),
                 [](const CSData& x) { return x.cmEnergy; });
  std::transform(measuredCS.begin(), measuredCS.end(),
                 _measuredCSData.cmEnergyError.data(),
                 [](const CSData& x) { return x.cmEnergyError; });
  std::transform(measuredCS.begin(), measuredCS.end(),
                 _measuredCSData.cs.data(),
                 [](const CSData& x) { return x.cs; });
  std::transform(measuredCS.begin(), measuredCS.end(),
                 _measuredCSData.csError.data(),
                 [](const CSData& x) { return x.csError; });
}

BaseISRSolver::~BaseISRSolver() {}

std::size_t BaseISRSolver::_getN() const { return _n; }

double BaseISRSolver::_energyThreshold() const { return _energyT; }

double BaseISRSolver::_sThreshold() const { return _energyT * _energyT; }

double BaseISRSolver::_s(std::size_t i) const {
  return _measuredCSData.cmEnergy(i) * _measuredCSData.cmEnergy(i);
}

double BaseISRSolver::_ecm(std::size_t i) const {
  return _measuredCSData.cmEnergy(i);
}

const Eigen::VectorXd& BaseISRSolver::_ecm() const {
  return _measuredCSData.cmEnergy;
}

Eigen::VectorXd& BaseISRSolver::_ecm() { return _measuredCSData.cmEnergy; }

const Eigen::VectorXd& BaseISRSolver::_ecmErr() const {
  return _measuredCSData.cmEnergyError;
}

Eigen::VectorXd& BaseISRSolver::_ecmErr() {
  return _measuredCSData.cmEnergyError;
}

const Eigen::VectorXd& BaseISRSolver::_vcs() const {
  return _measuredCSData.cs;
}

Eigen::VectorXd& BaseISRSolver::_vcs() { return _measuredCSData.cs; }

const Eigen::VectorXd& BaseISRSolver::_vcsErr() const {
  return _measuredCSData.csError;
}

Eigen::VectorXd& BaseISRSolver::_vcsErr() { return _measuredCSData.csError; }

Eigen::MatrixXd BaseISRSolver::_vcsInvErrMatrix() const {
  return _measuredCSData.csError.array().pow(-2.).matrix().asDiagonal();
}

const Eigen::VectorXd& BaseISRSolver::bcs() const { return _bornCS; }

Eigen::VectorXd& BaseISRSolver::_bcs() { return _bornCS; }
