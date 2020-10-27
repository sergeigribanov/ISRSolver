#ifndef __BASE_ISRSOLVER_HPP__
#define __BASE_ISRSOLVER_HPP__

#include <Eigen/Dense>
#include <string>
#include <TEfficiency.h>

#include "ISRSolverStructs.hpp"

typedef struct : std::exception {
  const char* what() const noexcept {
    return "[!] Wrong efficiency dimension.\n";
  }
} EfficiencyDimensionException;

class BaseISRSolver {
 public:
  BaseISRSolver(const std::string& inputPath, const InputOptions& inputOpts);
  BaseISRSolver(const BaseISRSolver&);
  virtual ~BaseISRSolver();
  virtual void solve() = 0;
  virtual void save(const std::string& outputPath,
                    const OutputOptions& outputOpts) = 0;
  const Eigen::VectorXd& bcs() const;
  const Eigen::VectorXd& ecm() const;
  const Eigen::VectorXd& ecmErr() const;
  const std::function<double(double, double)>& efficiency() const;

 protected:
  std::size_t _getN() const;
  double _energyThreshold() const;
  double _sThreshold() const;
  double _s(std::size_t) const;
  double _ecm(std::size_t) const;
  Eigen::VectorXd& _ecm();
  Eigen::VectorXd& _ecmErr();
  const Eigen::VectorXd& _vcs() const;
  Eigen::VectorXd& _vcs();
  const Eigen::VectorXd& _vcsErr() const;
  Eigen::VectorXd& _vcsErr();
  Eigen::VectorXd& _bcs();
  Eigen::MatrixXd _vcsInvErrMatrix() const; 
  void _setupEfficiency() noexcept(false);
  
 private:
  double _energyT;
  std::size_t _n;
  CSVecData _visibleCSData;
  std::function<double(double, double)> _efficiency;
  TEfficiency* _tefficiency;
  Eigen::VectorXd _bornCS;
};

#endif
