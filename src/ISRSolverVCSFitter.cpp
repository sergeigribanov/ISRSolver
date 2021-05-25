#include <algorithm>
#include "Integration.hpp"
#include "KuraevFadin.hpp"
#include "ISRSolverVCSFitter.hpp"

ISRSolverVCSFitFunction::ISRSolverVCSFitFunction(
    std::size_t n,
    double threshold,
    double* energy, double* vis_cs,
    double* energy_err, double* vis_cs_err,
    const std::function<double(double, const std::vector<double>&)>& fit_fcn,
    const std::function<double(double, double)>& eff_fcn) :
    _energySpread(false),
    _threshold(threshold),
    _errorDef(1.),
    _fcn(fit_fcn),
    _eff(eff_fcn) {
  _ecm.resize(n);
  _vcs.resize(n);
  if (energy_err) {
    _ecmErr.resize(n);
  } else {
    _ecmErr = std::vector<double>(n, 0.);
  }
  _vcsErr.resize(n);
  std::copy(energy, energy + n, _ecm.begin());
  std::copy(vis_cs, vis_cs + n, _vcs.begin());
  if (energy_err) {
    std::copy(energy_err, energy_err + n, _ecmErr.begin());
  }
  std::copy(vis_cs_err, vis_cs_err + n, _vcsErr.begin());
}

ISRSolverVCSFitFunction::~ISRSolverVCSFitFunction() {}

double ISRSolverVCSFitFunction::Up() const {
  return _errorDef;
}

double ISRSolverVCSFitFunction::operator()(
    const std::vector<double>& par) const {
  double chi2 = 0;
  std::function<double(double)> bcs_fcn =
      [par, this](double en) {
        const double result = this->_fcn(en, par);
        return result;
      };
  std::function<double(double)> vcs_no_spread =
      [bcs_fcn, this](double en) {
        const double sT = this->_threshold * this->_threshold;
        const double s = en * en;
        if (s <= sT) {
          return 0.;
        }
        const double result = convolutionKuraevFadin(
            en, bcs_fcn, 0, 1. - sT / s, this->_eff);
        return result;
      };
  for (std::size_t i = 0; i < _ecm.size(); ++i) {
    const double dvcs = gaussian_conv(
        _ecm[i], _ecmErr[i] * _ecmErr[i], vcs_no_spread) - _vcs[i];
    chi2 += dvcs * dvcs / _vcsErr[i] / _vcsErr[i];
  }
  return chi2;
}

void ISRSolverVCSFitFunction::setErrorDef(double def) {
  _errorDef = def;
}

void ISRSolverVCSFitFunction::enableEnergySpread() {
  _energySpread = true;
}

void ISRSolverVCSFitFunction::disableEnergySpread() {
  _energySpread = false;
}
