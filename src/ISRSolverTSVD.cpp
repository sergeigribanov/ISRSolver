#include <iostream>
#include <Eigen/Core>
#include <Eigen/SVD>
#include "ISRSolverTSVD.hpp"

ISRSolverTSVD::ISRSolverTSVD(const std::string& inputPath,
                             const InputOptions& inputOpts,
                             int truncIndexUpperLimit):
    ISRSolverSLAE(inputPath, inputOpts),
    _truncIndexUpperLimit(truncIndexUpperLimit) {}

ISRSolverTSVD::ISRSolverTSVD(const ISRSolverTSVD& solver):
    ISRSolverSLAE::ISRSolverSLAE(solver),
    _truncIndexUpperLimit(solver._truncIndexUpperLimit),
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
  Eigen::MatrixXd mK = _mU.block(0, 0, _mU.rows(), _truncIndexUpperLimit) *
                       _mSing.head(_truncIndexUpperLimit).asDiagonal() *
                       _mV.block(0, 0, _mV.rows(), _truncIndexUpperLimit).transpose();
  _bcs() = mK.completeOrthogonalDecomposition().solve(_vcs());
  _getInverseBornCSErrorMatrix() = mK.transpose() * _vcsInvErrMatrix() * mK;
}

void ISRSolverTSVD::setTruncIndexUpperLimit(int truncIndexUpperLimit) {
  _truncIndexUpperLimit = truncIndexUpperLimit;
}

int ISRSolverTSVD::getTruncIndexUpperLimit() const {
  return _truncIndexUpperLimit;
}
