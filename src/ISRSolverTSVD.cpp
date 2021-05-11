#include <iostream>
#include <Eigen/Core>
#include <Eigen/SVD>
#include "ISRSolverTSVD.hpp"

ISRSolverTSVD::ISRSolverTSVD(
    std::size_t numberOfPoints,
    double* energy, double* visibleCS,
    double* energyErr, double* visibleCSErr,
    double thresholdEnergy,
    const std::function<double(double, double)>&
    efficiency) :
    ISRSolverSLE(numberOfPoints,
                 energy, visibleCS,
                 energyErr, visibleCSErr,
                 thresholdEnergy,
                 efficiency),
    _upperTSVDIndex(numberOfPoints),
    _keepOne(false) {}

ISRSolverTSVD::ISRSolverTSVD(TGraphErrors* vcsGraph,
                             double thresholdEnergy,
                             int upperTSVDIndex) :
    ISRSolverSLE(vcsGraph, thresholdEnergy),
    _upperTSVDIndex(upperTSVDIndex),
    _keepOne(false) {}

ISRSolverTSVD::ISRSolverTSVD(TGraphErrors* vcsGraph,
                             TEfficiency* eff,
                             double thresholdEnergy,
                             int upperTSVDIndex) :
    ISRSolverSLE(vcsGraph, eff, thresholdEnergy),
    _upperTSVDIndex(upperTSVDIndex),
    _keepOne(false) {}

ISRSolverTSVD::ISRSolverTSVD(const std::string& inputPath,
                             const InputOptions& inputOpts,
                             int upperTSVDIndex) :
    ISRSolverSLE(inputPath, inputOpts),
    _upperTSVDIndex(upperTSVDIndex),
    _keepOne(false) {}

ISRSolverTSVD::ISRSolverTSVD(const ISRSolverTSVD& solver):
    ISRSolverSLE::ISRSolverSLE(solver),
    _upperTSVDIndex(solver._upperTSVDIndex),
    _keepOne(solver._keepOne),
    _mU(solver._mU),
    _mV(solver._mV),
    _mSing(solver._mSing) {}

ISRSolverTSVD::~ISRSolverTSVD() {}

void ISRSolverTSVD::solve() {
  if (!_isEqMatrixPrepared) {
    evalEqMatrix();
    _isEqMatrixPrepared = true;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(getIntegralOperatorMatrix(),
                                          Eigen::ComputeFullV | Eigen::ComputeFullU);
    _mU = svd.matrixU();
    _mV = svd.matrixV();
    _mSing = svd.singularValues();
  }
  int firstIndex = 0;
  int n = _upperTSVDIndex;
  if (_keepOne) {
    firstIndex = _upperTSVDIndex - 1;
    n = 1;
  }
  Eigen::MatrixXd mK = _mU.block(0, firstIndex, _mU.rows(), n) *
                       _mSing.segment(firstIndex, n).asDiagonal() *
                       _mV.block(0, firstIndex, _mV.rows(), n).transpose();
  _bcs() = mK.completeOrthogonalDecomposition().solve(_vcs());
  _getBornCSCovMatrix() = (mK.transpose() * _vcsInvCovMatrix() * mK).inverse();
}

void ISRSolverTSVD::setUpperTSVDIndex(int upperTSVDIndex) {
  _upperTSVDIndex = upperTSVDIndex;
}

int ISRSolverTSVD::getUpperTSVDIndex() const {
  return _upperTSVDIndex;
}

void ISRSolverTSVD::enableKeepOne() {
  _keepOne = true;
}

void ISRSolverTSVD::disableKeepOne() {
  _keepOne = false;
}
