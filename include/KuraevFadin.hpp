#ifndef _KURAEV_FADIN_HPP_
#define _KURAEV_FADIN_HPP_
#include <functional>

/**
 * Convolution of a function with the Kuraev-Fadin kernel function and a detection efficiency
 * @param energy a center-of-mass energy
 * @param fcn a function that is convoluted
 * @param min_x a lower integration limit
 * @param max_x an upper integration limut
 * @param efficiency a detection efficiency (default value = 1)
 */
double convolutionKuraevFadin(double energy,
                              const std::function<double(double)>& fcn,
                              double min_x,
                              double max_x,
                              const std::function<double(double, double)>& efficiency =
                              [](double, double) {return 1.;});

double convolutionKuraevFadin_1d(double energy,
                              const std::function<double(double)>& fcn,
                              double min_x,
                              double max_x,
                              const std::function<double(double)>& efficiency =
                              [](double) {return 1.;});

/**
 * The Kuraev-Fadin kernel function.
 * @param x an argument x
 * @param s a square of a center-of-mass energy
 */
double kernelKuraevFadin(double x, double s);

#endif
