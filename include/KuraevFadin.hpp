#pragma once

#include <functional>


/**
 * @brief Convolve function with Kuraev-Fadin kernel (2D efficiency)
 * @param energy Center-of-mass energy
 * @param fcn Function to convolute
 * @param min_x Lower integration limit
 * @param max_x Upper integration limit
 * @param efficiency Detection efficiency (default=1)
 */
double convolutionKuraevFadin(
    double energy,
    const std::function<double(double)>& fcn,
    double min_x,
    double max_x,
    const std::function<double(double, double)>& efficiency = 
        [](double, double) { return 1.; });

/**
 * @brief Convolve function with Kuraev-Fadin kernel (1D efficiency)
 * @param energy Center-of-mass energy
 * @param fcn Function to convolute
 * @param min_x Lower integration limit
 * @param max_x Upper integration limit
 * @param efficiency Detection efficiency (default=1)
 */
double convolutionKuraevFadin_1d(
    double energy,
    const std::function<double(double)>& fcn,
    double min_x,
    double max_x,
    const std::function<double(double)>& efficiency = 
        [](double) { return 1.; });

/**
 * @brief Kuraev-Fadin kernel function
 * @param x Energy fraction
 * @param s Squared center-of-mass energy
 */
double kernelKuraevFadin(double x, double s);

/**
 * @brief Beta function for ISR calculations
 * @param s Squared center-of-mass energy
 */
double fBeta(double s);

/**
 * @brief D-function for ISR calculations
 * @param mX 1 - x (x = photon energy fraction)
 * @param s Squared center-of-mass energy
 */
double funcD(double mX, double s);

/**
 * @brief Integral of F kernel over x-range
 * @param x_min Lower integration limit
 * @param x_max Upper integration limit
 * @param s Squared center-of-mass energy
 */
double funcFIntegral(double x_min, double x_max, double s);

/**
 * @brief Majorant function
 * @param mX 1 - x
 * @param s Squared center-of-mass energy
 * @param koeff Scaling factor
 */
double majorantDFunction(double mX, double s, double koeff);
