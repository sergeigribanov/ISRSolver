#ifndef _ISRSOLVER_VCS_FITTER_HPP_
#define _ISRSOLVER_VCS_FITTER_HPP_
#include <Minuit2/FCNBase.h>
#include <Minuit2/FunctionMinimizer.h>
#include <Minuit2/FunctionMinimum.h>
#include <functional>
#include <string>
#include <vector>

class ISRSolverVCSFitFunction : public ROOT::Minuit2::FCNBase {
 public:
  ISRSolverVCSFitFunction(std::size_t n,
                          double threshold,
                          double* energy, double* vis_cs,
                          double* energy_err, double* vis_cs_err,
                          const std::function<double(double, const std::vector<double>&)>& fit_fcn,
                          const std::function<double(double, double)>& eff_fcn =
                          [](double, double) {return 1.;});
  virtual ~ISRSolverVCSFitFunction();
  virtual double Up() const override final;
  virtual double operator()(const std::vector<double>&) const override final;
  void setErrorDef(double def);
  void enableEnergySpread();
  void disableEnergySpread();
  void saveResults(const ROOT::Minuit2::FunctionMinimum&, const std::string&) const;
 private:
  bool _energySpread;
  double _threshold;
  double _errorDef;
  std::function<double(double, const std::vector<double>&)> _fcn;
  std::function<double(double, double)> _eff;
  std::vector<double> _ecm;
  std::vector<double> _ecmErr;
  std::vector<double> _vcs;
  std::vector<double> _vcsErr;
};

#endif
