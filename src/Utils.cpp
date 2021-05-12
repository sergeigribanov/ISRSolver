#include <nlopt.hpp>
#include <TRandom3.h>
#include "ISRSolverTikhonov.hpp"
#include "Utils.hpp"

Eigen::VectorXd randomDrawVisCS(const Eigen::VectorXd& vcs,
                                const Eigen::VectorXd& vcsErr) {
  Eigen::VectorXd result = vcs;
  for (int i = 0; i < result.size(); ++i) {
    result(i) = gRandom->Gaus(result(i), vcsErr(i));
  }
  return result;
}

/**
 * Objective function that return the L-curve curvature with a negative sign
 */
double lambdaObjective(unsigned n, const double* plambda, double* grad, void* solver) {
  auto sp = reinterpret_cast<ISRSolverTikhonov*>(solver);
  /**
   * Setting regularization parameter
   */
  sp->setLambda(*plambda);
  /**
   * Finding a numerical solution
   */
  sp->solve();
  if (grad) {
    /**
     * Evaluate L-curve curvature gradient
     */
    grad[0] = sp->evalLCurveCurvatureDerivative();
  }
  /**
   * Return L-curve curvature
   */
  return sp->evalLCurveCurvature();
}
