#ifndef _INTEGRATION_HPP_
#define _INTEGRATION_HPP_
#include <functional>

/**
 * Adaptive integration using GSL
 * @param fcn an integrand
 * @parm a a lower integration limit
 * @param b an upper integration limit
 * @param error an integration absolute error
 */
double integrate(std::function<double(double)>& fcn, double a, double b, double& error);
/**
 * Adaptive singular integration using GSL
 * @param fcn an integrand
 * @param a a lower integration limit
 * @param b an upper integration limit
 * @param error an integration absolute error
 */
double integrateS(std::function<double(double)>& fcn, double a, double b, double& error);
/**
 * Gaussian convolution
 * @param energy a mean center-of-mass energy
 * @param sigma2 a square of standard deviation for center-of-mass energy
 * @param fcn an integrand
 */
double gaussian_conv(double energy, double sigma2, std::function<double(double)>& fcn);

#endif
