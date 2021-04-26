#ifndef _KURAEV_FADIN_HPP_
#define _KURAEV_FADIN_HPP_
#include <functional>

double kuraev_fadin_convolution(double, const std::function<double(double)>&,
                                double, double,
				const std::function<double(double, double)>& =
				[](double, double) {return 1.;});

double kuraev_fadin_kernel(double, double);

#endif
