#ifndef __GET_INTEGRAL_H__
#define __GET_INTEGRAL_H__
#include <functional>

double getIntegral(std::function<double(double)>&, double, double, double&);
double getIntegralS(std::function<double(double)>&, double, double, double&);

#endif
