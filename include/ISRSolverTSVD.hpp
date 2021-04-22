#ifndef __ISRSOLVER_TSVD_HPP__
#define __ISRSOLVER_TSVD_HPP__

#include "ISRSolverSLE.hpp"

class ISRSolverTSVD : public ISRSolverSLE {
 public:
  ISRSolverTSVD(const std::string&, const InputOptions&, int);
  ISRSolverTSVD(const ISRSolverTSVD&);
  virtual ~ISRSolverTSVD();
  virtual void solve() override;
  void setTruncIndexUpperLimit(int);
  int getTruncIndexUpperLimit() const;
  void enableKeepOne();
  void disableKeepOne();
 private:
  int _truncIndexUpperLimit;
  bool _keepOne;
  Eigen::MatrixXd _mU;
  Eigen::MatrixXd _mV;
  Eigen::VectorXd _mSing;
};

#endif
