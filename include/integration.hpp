#ifndef __INTEGRATION_HPP__
#define __INTEGRATION_HPP__
#include <functional>

double getIntegral(std::function<double(double)>&, double, double, double&);
double getIntegralS(std::function<double(double)>&, double, double, double&);

#endif
