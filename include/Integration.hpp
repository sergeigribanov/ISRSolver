#ifndef _INTEGRATION_HPP_
#define _INTEGRATION_HPP_
#include <functional>

double integrate(std::function<double(double)>&, double, double, double&);
double integrateS(std::function<double(double)>&, double, double, double&);
double gaussian_conv(double, double, std::function<double(double)>&);

#endif
