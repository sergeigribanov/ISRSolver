#include <TGraphErrors.h>
#include <algorithm>
#include "Integration.hpp"
#include "KuraevFadin.hpp"
#include "ISRSolverVCSFitter.hpp"

#include <iostream>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>

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
  std::vector<std::size_t> inds(n);
  std::iota(inds.begin(), inds.end(), 0);
  std::sort(inds.begin(), inds.end(),
            [this](std::size_t i1, std::size_t i2) {
              return (this->_ecm)[i1] < (this->_ecm)[i2];});
  auto sort_v = [n, inds](std::vector<double> &v) {
    std::vector<double> tmp_v(n);
    std::transform(inds.begin(), inds.end(), tmp_v.begin(),
                   [v](std::size_t i) { return v[i]; });
    v = tmp_v;
  };
  sort_v(_ecm);
  sort_v(_vcs);
  sort_v(_ecmErr);
  sort_v(_vcsErr);
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
            en, bcs_fcn, 0., 1. - sT / s, this->_eff);
        return result;
      };
  for (std::size_t i = 0; i < _ecm.size(); ++i) {
    double dvcs = 0.;
    if (_energySpread) {
      dvcs = gaussian_conv(_ecm[i], _ecmErr[i] * _ecmErr[i], vcs_no_spread) - _vcs[i];
    } else {
      dvcs = vcs_no_spread(_ecm[i]) - _vcs[i];
    }
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

void ISRSolverVCSFitFunction::saveResults(
    const ROOT::Minuit2::FunctionMinimum &fcn_min,
    const std::string &path) const {
  auto params = fcn_min.UserParameters().Params();
  std::function<double(double)> bcs = [params, this](double en) {
    const double result = this->_fcn(en, params);
    return result;
  };
  std::function<double(double *, double *)> bcs_fcn = [bcs](double* x, double*) {
    const double result = bcs(x[0]);
    return result;
  };
  std::function<double(double)> vcs_no_spread =
      [bcs, this](double en) {
        const double sT = this->_threshold * this->_threshold;
        const double s = en * en;
        if (s <= sT) {
          return 0.;
        }
        const double result =
            convolutionKuraevFadin(en, bcs, 0., 1. - sT / s, this->_eff);
        return result;
      };

  std::function<double(double*, double*)> vcs_no_spread_fcn = [vcs_no_spread](double* x, double*) {
    const double result = vcs_no_spread(x[0]);
    return result;
  };
  TF1 f_bcs("f_bcs", bcs_fcn, *std::min_element(_ecm.begin(), _ecm.end()),
            *std::max_element(_ecm.begin(), _ecm.end()),  0);
  f_bcs.SetNpx(100000);
  TF1 f_vcs("f_vcs", vcs_no_spread_fcn, *std::min_element(_ecm.begin(), _ecm.end()),
            *std::max_element(_ecm.begin(), _ecm.end()),  0);
  f_vcs.SetNpx(100000);
  auto fl = TFile::Open(path.c_str(), "recreate");
  fl->cd();
  f_bcs.Write();
  f_vcs.Write();
  std::vector<double> vcs_gc;
  if (_energySpread) {
    for (std::size_t i = 0; i < _ecm.size(); ++i) {
      vcs_gc.push_back(gaussian_conv(_ecm[i], _ecmErr[i] * _ecmErr[i], vcs_no_spread));
    }
  } else {
    for (std::size_t i = 0; i <  _ecm.size(); ++i) {
      vcs_gc.push_back(f_vcs.Eval(_ecm[i]));
    }
  }
  TGraph gr_vcs_gauss_conv(_ecm.size(), _ecm.data(), vcs_gc.data());
  gr_vcs_gauss_conv.SetMarkerStyle(20);
  gr_vcs_gauss_conv.Write("vcs_gauss_conv");
  std::vector<double> rad_delta;
  std::vector<double> vec_bcs;
  std::vector<double> vec_bcs_err;
  for (std::size_t i = 0; i < _ecm.size(); ++i) {
    double val = vcs_gc[i] / f_bcs.Eval(_ecm[i]) - 1.;
    rad_delta.push_back(val);
    vec_bcs.push_back(_vcs[i] / (val + 1.));
    vec_bcs_err.push_back(_vcsErr[i] / (val + 1.));
  }
  TGraph gr_rad_delta(_ecm.size(), _ecm.data(), rad_delta.data());
  gr_rad_delta.SetMarkerStyle(20);
  gr_rad_delta.Write("rad_delta");
  TGraphErrors gr_vcs(_ecm.size(), _ecm.data(), _vcs.data(), _ecmErr.data(), _vcsErr.data());
  gr_vcs.SetMarkerStyle(20);
  gr_vcs.Write("vcs");
  TGraphErrors gr_bcs(_ecm.size(), _ecm.data(), vec_bcs.data(), 0, vec_bcs_err.data());
  gr_bcs.Write("bcs");
  fl->Close();
  delete fl;
}
