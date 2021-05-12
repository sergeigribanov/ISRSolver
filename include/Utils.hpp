#ifndef _UTILS_HPP_
#define _UTILS_HPP_
#include <iostream>
#include <string>
#include <Eigen/Dense>
#include <TFile.h>

/**
 * Random redraw of a visible cross section
 * @param vcs an initial visible cross section
 * @param vcsErr a visible cross section errors
 */
Eigen::VectorXd randomDrawVisCS(const Eigen::VectorXd& vcs,
                                const Eigen::VectorXd& vcsErr);

double lambdaObjective(unsigned n, const double* plambda, double* grad, void* solver);

#endif
