#include <functional>
#include <nlohmann/json.hpp>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMatrixD.h>
#include "IterISRInterpSolver.hpp"
using json = nlohmann::json;

IterISRInterpSolver::IterISRInterpSolver(const std::string& inputPath,
                                         const InputOptions& inputOpts) :
    ISRSolverSLE(inputPath, inputOpts),
    _nIter(10) {}

IterISRInterpSolver::IterISRInterpSolver(const IterISRInterpSolver& solver) :
    ISRSolverSLE(solver),
    _nIter(solver._nIter),
    _radcorr(solver._radcorr) {}

IterISRInterpSolver::~IterISRInterpSolver() {}

void IterISRInterpSolver::solve() {
  if (!_isEqMatrixPrepared) {
    evalEqMatrix();
    _isEqMatrixPrepared = true;
  }
  _bcs() = _vcs();
  for (std::size_t iter = 0; iter < _nIter; ++iter) {
    Eigen::VectorXd curVCS = getIntegralOperatorMatrix() * _bcs();
    _radcorr = curVCS.array() / _bcs().array();
    _bcs() = _vcs().array() / _radcorr.array();
    _getBornCSCovMatrix() = (_vcsErr().array() / _radcorr.array()).pow(2.).matrix().asDiagonal();
  }
}

void IterISRInterpSolver::setNumOfIters(std::size_t nIter) {
  _nIter = nIter;
}

std::size_t IterISRInterpSolver::getNumOfIters() const {
  return _nIter;
}

void IterISRInterpSolver::save(const std::string& outputPath,
                               const OutputOptions& outputOpts) {
  TGraphErrors gvcs(_getN(), _ecm().data(), _vcs().data(), _ecmErr().data(),
                   _vcsErr().data());
  Eigen::VectorXd bcsErr = getBornCSCovMatrix().diagonal().array().pow(0.5);
  TGraphErrors gbcs(_getN(), _ecm().data(), _bcs().data(), 0, bcsErr.data());
  TGraph gradcor(_getN(), _ecm().data(), _radcorr.data());
  TMatrixD intergalOperatorMatrix(_getN(), _getN());
  Eigen::MatrixXd tmpIntOpM = getIntegralOperatorMatrix().transpose();
  intergalOperatorMatrix.SetMatrixArray(tmpIntOpM.data());
  TMatrixD bornCSCovMatrix(_getN(), _getN());
  Eigen::MatrixXd tmpCovM = getBornCSCovMatrix().transpose();
  bornCSCovMatrix.SetMatrixArray(tmpCovM.data());
  TMatrixD bornCSInvCovMatrix(_getN(), _getN());
  Eigen::MatrixXd tmpInvCovM = getBornCSCovMatrix().inverse().transpose();
  bornCSInvCovMatrix.SetMatrixArray(tmpInvCovM.data());
  auto fl = TFile::Open(outputPath.c_str(), "recreate");
  fl->cd();
  gvcs.Write(outputOpts.visibleCSGraphName.c_str());
  gbcs.Write(outputOpts.bornCSGraphName.c_str());
  gradcor.Write("radiative_correction");
  intergalOperatorMatrix.Write("intergalOperatorMatrix");
  bornCSCovMatrix.Write("covMatrixBornCS");
  bornCSInvCovMatrix.Write("invCovMatrixBornCS");
  fl->Close();
  delete fl;
}
