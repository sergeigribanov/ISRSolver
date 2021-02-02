#include <iostream>
#include <Eigen/Core>
#include <Eigen/SVD>
#include "ISRSolverTSVD.hpp"

ISRSolverTSVD::ISRSolverTSVD(const std::string& inputPath,
                             const InputOptions& inputOpts,
                             int truncIndexUpperLimit):
    ISRSolverSLAE(inputPath, inputOpts),
    _truncIndexUpperLimit(truncIndexUpperLimit),
    _keepOne(false) {}

ISRSolverTSVD::ISRSolverTSVD(const ISRSolverTSVD& solver):
    ISRSolverSLAE::ISRSolverSLAE(solver),
    _truncIndexUpperLimit(solver._truncIndexUpperLimit),
    _keepOne(solver._keepOne),
    _mU(solver._mU),
    _mV(solver._mV),
    _mSing(solver._mSing) {}

ISRSolverTSVD::~ISRSolverTSVD() {}

void ISRSolverTSVD::solve() {
  if (!_isEqMatrixPrepared) {
    _evalEqMatrix();
    _isEqMatrixPrepared = true;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(getIntegralOperatorMatrix(),
                                          Eigen::ComputeFullV | Eigen::ComputeFullU);
    _mU = svd.matrixU();
    _mV = svd.matrixV();
    _mSing = svd.singularValues();
  }
  int firstIndex = 0;
  int n = _truncIndexUpperLimit;
  if (_keepOne) {
    firstIndex = _truncIndexUpperLimit - 1;
    n = 1;
  }
  Eigen::MatrixXd mK = _mU.block(0, firstIndex, _mU.rows(), n) *
                       _mSing.segment(firstIndex, n).asDiagonal() *
                       _mV.block(0, firstIndex, _mV.rows(), n).transpose();
  _bcs() = mK.completeOrthogonalDecomposition().solve(_vcs());
  _getInverseBornCSErrorMatrix() = mK.transpose() * _vcsInvErrMatrix() * mK;
}

void ISRSolverTSVD::setTruncIndexUpperLimit(int truncIndexUpperLimit) {
  _truncIndexUpperLimit = truncIndexUpperLimit;
}

int ISRSolverTSVD::getTruncIndexUpperLimit() const {
  return _truncIndexUpperLimit;
}

void ISRSolverTSVD::enableKeepOne() {
  _keepOne = true;
}

void ISRSolverTSVD::disableKeepOne() {
  _keepOne = false;
}
