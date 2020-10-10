#ifndef __BASE_ISRSOLVER_HPP__
#define __BASE_ISRSOLVER_HPP__

#include <Eigen/Dense>
#include <string>

#include "ISRSolverStructs.hpp"

class BaseISRSolver {
 public:
  BaseISRSolver(const std::string& inputPath, const InputOptions& inputOpts);
  virtual ~BaseISRSolver();
  virtual void solve() = 0;
  virtual void save(const std::string& outputPath,
                    const OutputOptions& outputOpts) = 0;
  const Eigen::VectorXd& bcs() const;
  const Eigen::VectorXd& ecm() const;

 protected:
  std::size_t _getN() const;
  double _energyThreshold() const;
  double _sThreshold() const;
  double _s(std::size_t) const;
  double _ecm(std::size_t) const;
  const Eigen::VectorXd& _ecm() const;
  Eigen::VectorXd& _ecm();
  const Eigen::VectorXd& _ecmErr() const;
  Eigen::VectorXd& _ecmErr();
  const Eigen::VectorXd& _vcs() const;
  Eigen::VectorXd& _vcs();
  const Eigen::VectorXd& _vcsErr() const;
  Eigen::VectorXd& _vcsErr();
  Eigen::VectorXd& _bcs();
  Eigen::MatrixXd _vcsInvErrMatrix() const;

 private:
  double _energyT;
  std::size_t _n;
  CSVecData _visibleCSData;
  Eigen::VectorXd _bornCS;
};

#endif
