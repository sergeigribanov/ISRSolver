#ifndef __KURAEV_FADIN_HPP__
#define __KURAEV_FADIN_HPP__
#include <TGraphErrors.h>

#include <functional>

double kuraev_fadin_polinomial_convolution(double, double, double, int);

double kuraev_fadin_convolution(double, const std::function<double(double)>&,
                                double, double);

double kuraev_fadin_kernel(double, double);

#endif
