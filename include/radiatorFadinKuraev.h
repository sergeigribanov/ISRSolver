#ifndef __RADIATOR_FADIN_KURAEV_H__
#define __RADIATOR_FADIN_KURAEV_H__
#include "TGraphErrors.h"
#include <functional>

double radIntegral(double, double, double, int);

double getFadinIntegral(double,
                        const std::function<double(double)>&,
                        double);

double radiatorFadinKuraev(double, double);

#endif
