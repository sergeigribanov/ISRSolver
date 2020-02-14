#ifndef __KURAEV_FADIN_HPP__
#define __KURAEV_FADIN_HPP__
#include <TGraphErrors.h>

#include <functional>

double radIntegral(double, double, double, int);

double getFadinIntegral(double, const std::function<double(double)>&, double,
                        double);

double radiatorFadinKuraev(double, double);

#endif
