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
 * @see Chi2TestArgs
 */
void chi2Test(ISRSolverSLE* solver,
              const Chi2TestArgs& args,
              const Eigen::VectorXd& vcs,
              const Eigen::VectorXd& bcs0,
              const Eigen::VectorXd& vcsErr);

/**
 * Chi-square model test
 *
 */
void chi2TestModel(ISRSolverSLE*,
                   const Chi2TestModelArgs&);

void chi2TestData(ISRSolverSLE*,
                  const Chi2TestArgs&);

#endif
