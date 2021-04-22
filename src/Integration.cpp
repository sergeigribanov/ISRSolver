#define _USE_MATH_DEFINES
#include "Integration.hpp"
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include <iostream>

double wrapper(double x, void* fcnp) {
  auto fp = static_cast<std::function<double(double)>*>(fcnp);
  return (*fp)(x);
}

double integrateS(std::function<double(double)>& fcn, double a, double b,
                  double& error) {
  int N = 100000;
  gsl_error_handler_t* old_handler = gsl_set_error_handler_off();
  gsl_integration_workspace* w = gsl_integration_workspace_alloc(N);
  gsl_function F;
  F.function = &wrapper;
  F.params = &fcn;
  double result;
  double relerr = 1.0e-12;
  int status = 1;
  while (status) {
    status = gsl_integration_qags(&F, a, b, 1.e-12, relerr, N, w, &result, &error);

    if (relerr < 1.e-3) {
      relerr *= 10;
    } else {
      relerr *= 1.1;
      if (status) {
        std::cout << "AS: Warning: tolerance increased to " << relerr
                  << std::endl;
      }
    }
  }

  gsl_integration_workspace_free(w);
  gsl_set_error_handler(old_handler);

  return result;
}

double integrate(std::function<double(double)>& fcn, double a, double b,
                 double& error) {
  int N = 1000000;
  gsl_error_handler_t* old_handler = gsl_set_error_handler_off();
  gsl_integration_workspace* w = gsl_integration_workspace_alloc(N);
  gsl_function F;
  F.function = &wrapper;
  F.params = &fcn;
  double result;
  double relerr = 1.0e-12;
  int status = 1;
  while (status) {
    status =
        gsl_integration_qag(&F, a, b, 1.e-12, relerr, N, 6, w, &result, &error);

    if (relerr < 1.e-3) {
      relerr *= 10;
    } else {
      relerr *= 1.1;
      if (status) {
        std::cout << "A: Warning: tolerance increased to " << relerr
                  << std::endl;
      }
    }
  }
  gsl_integration_workspace_free(w);
  gsl_set_error_handler(old_handler);
  return result;
}

double gaussian_conv(double energy,
                     double sigma2,
                     std::function<double(double)>& fcn) {
  const double b = 0.5 / sigma2;
  gsl_integration_fixed_workspace * w;
  const gsl_integration_fixed_type * T = gsl_integration_fixed_hermite;
  w = gsl_integration_fixed_alloc(T, 6, energy, b, 0., 0.);
  gsl_function F;
  F.function = &wrapper;
  F.params = &fcn;
  double result;
  gsl_integration_fixed(&F, &result, w);
  result /= std::sqrt(2 * M_PI * sigma2);
  return result;
}
