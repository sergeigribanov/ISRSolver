#ifndef _CHI2_TEST_HPP_
#define _CHI2_TEST_HPP_
#include <string>
#include "ISRSolverStructs.hpp"
#include "ISRSolverSLE.hpp"

/**
 * Chi-square test
 * @param solver a solver
 * @param args a chi-square test arguments
 * @param vcs a visible cross section
 * @param bcs0 a reference Born cross section
 * @param vcsErr a visible cross section errors
 * @see ISRSolverSLE
 * @see Chi2TestArgs
 */
void chi2Test(ISRSolverSLE* solver,
              const Chi2TestArgs& args,
              const Eigen::VectorXd& vcs,
              const Eigen::VectorXd& bcs0,
              const Eigen::VectorXd& vcsErr);

/**
 * Chi-square model test
 * @parm solver a solver
 * @parm modelArgs a structure with input arguments
 * @see ISRSolverSLE
 * @see Chi2TestModelArgs
 */
void chi2TestModel(ISRSolverSLE* solver,
                   const Chi2TestModelArgs& modelArgs);

/**
 * Chi-square data test
 * @param solver a solver
 * @param args a structure with input arguments
 * @see ISRSolverSLE
 * @see Chi2TestArgs
 */
void chi2TestData(ISRSolverSLE* solver,
                  const Chi2TestArgs& args);

#endif
