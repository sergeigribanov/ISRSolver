#ifndef __INTEGRATION_HPP__
#define __INTEGRATION_HPP__
#include <functional>

double integrate(std::function<double(double)>&, double, double, double&);
double integrateS(std::function<double(double)>&, double, double, double&);

#endif
