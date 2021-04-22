#ifndef _ITER_ISR_INTERP_SOLVER_HPP_
#define _ITER_ISR_INTERP_SOLVER_HPP_
#include "ISRSolverSLE.hpp"
#include "Interpolator.hpp"

class IterISRInterpSolver : public ISRSolverSLE {
 public:
  IterISRInterpSolver(const std::string&, const InputOptions&);
  IterISRInterpSolver(const IterISRInterpSolver&);
  virtual ~IterISRInterpSolver();
  virtual void solve() override;
  void setNumOfIters(std::size_t);
  std::size_t getNumOfIters() const;
  virtual void save(const std::string&, const OutputOptions&) override;
 private:
  std::size_t _nIter;
  Eigen::VectorXd _radcorr;
};

#endif
