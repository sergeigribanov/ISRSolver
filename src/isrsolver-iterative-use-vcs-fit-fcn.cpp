#include <iostream>
#include <string>
#include <vector>
#include <functional>
#include <limits>
#include <boost/program_options.hpp>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <cmath>
#include <Eigen/Dense>
#include <TEfficiency.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TFile.h>
#include "Integration.hpp"
#include "KuraevFadin.hpp"
namespace po = boost::program_options;

/**
 * A part of program options
 */
typedef struct {
  /**
   * Number of iterations
   */
  std::size_t niter;
  /**
   * Number of points at which the radiative correction is calculated
   */
  std::size_t n;
  /**
   * Threshold energy
   */
  double thsd;
  /**
   * Center-of-mass energy spread
   */
  double energy_spread;
  /**
   * Name of the visible cross section graph (TGraphErrors)
   */
  std::string vcs_name;
  /**
   * Name of the visible cross section fit function (TF1)
   */
  std::string fcn_vcs_name;
  /**
   * Name of the detection efficiency object (TEfficiency)
   */
  std::string efficiency_name;
  /**
   * Path to the input file that contains the visible cross section
   * and the detection efficiency
   */
  std::string ifname;
  /**
   * Path to the output file
   */
  std::string ofname;
} CmdOptions;

/**
 * Setting up program options
 */
void setOptions(po::options_description* desc, CmdOptions* opts) {
  desc->add_options()
      ("help,h", "help message")
      ("enable-energy-spread,g", "enable energy spread")
      ("niter,n", po::value<std::size_t>(&(opts->niter))->default_value(10),
       "number of iterations")
      ("number-of-points,p", po::value<std::size_t>(&(opts->n))->default_value(100),
       "number of iterations")
      ("thsd,t", po::value<double>(&(opts->thsd)), "threshold (GeV)")
      ("energy-spread,s", po::value<double>(&(opts->energy_spread))->default_value(0.),
       "Center-of-mass energy spread (GeV)")
      ("vcs-name,v", po::value<std::string>(&(opts->vcs_name))->default_value("vcs"),
       "name of the visible cross section graph (TGraphErrors*)")
      ("fcn-vcs-name,f", po::value<std::string>(&(opts->fcn_vcs_name))->default_value("fcn_vcs"),
       "name of the visible cross section fit function (TF1*)")
      ("efficiency-name,e", po::value<std::string>(&(opts->efficiency_name)),
       "name of the detection efficiency object (TEfficiency*)")
      ( "ifname,i",
        po::value<std::string>(&(opts->ifname))->default_value("vcs.root"),
        "Path to input file.")
      ("ofname,o",
       po::value<std::string>(&(opts->ofname))->default_value("bcs.root"),
       "Path to output file.");
}

/**
 * Help message
 */
void help(const po::options_description& desc) {
  std::cout << desc << std::endl;
}

int main(int argc, char* argv[]) {
  po::options_description desc("   Common used iterative solver (Using fit function for the visible cross section "
                               "and same energy spread at each center-of-mass energy point). Allowed options");
  CmdOptions opts;
  setOptions(&desc, &opts);
  po::variables_map vmap;
  po::store(po::parse_command_line(argc, argv, desc), vmap);
  po::notify(vmap);
  if (vmap.count("help")) {
    help(desc);
    return 0;
  }
  TEfficiency* teff = nullptr;
  /**
   * Opening the input file
   */
  auto fl = TFile::Open(opts.ifname.c_str(), "read");
  /**
   * Reading the visible cross section
   */
  auto vcs = dynamic_cast<TGraphErrors*>(fl->Get(opts.vcs_name.c_str())->Clone());
  /**
   * Reading the visible cross section fit function
   */
  auto fcn = dynamic_cast<TF1*>(fl->Get(opts.fcn_vcs_name.c_str())->Clone());
  /**
   * Reading the detection efficiency
   */
  if (vmap.count("efficiency-name")) {
    teff = dynamic_cast<TEfficiency*>(fl->Get(opts.efficiency_name.c_str())->Clone());
  }
  fl->Close();
  delete fl;
  Eigen::VectorXd ecm = Eigen::VectorXd::Zero(opts.n);
  const double minen = fcn->GetXmin();
  const double maxen = fcn->GetXmax();
  const double eh = (maxen - opts.thsd) / (opts.n - 1);
  for (std::size_t i = 0; i < opts.n; ++i) {
    ecm(i) = opts.thsd + eh * i;
  }
  Eigen::VectorXd radCorr = Eigen::VectorXd::Zero(opts.n);
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_linear, opts.n);
  gsl_spline_init (spline, ecm.data(), radCorr.data(), opts.n);
  std::function<double(double)> radFCN =
      [&spline, &acc](double en) {
        double result = gsl_spline_eval(spline, en, acc);
        return result;
      };
  /**
   * Converting the function to a form of std::function
   */
  std::function<double(double)> born_fcn =
      [opts, maxen, &fcn, &radFCN](double en) {
    const double e0 = opts.thsd;
    const double e1 = fcn->GetXmin();
    if (en <= e0) {
      return 0.;
    }
    double result = 0;
    if (en >= maxen) {
      result = fcn->Eval(maxen) / (1. + radFCN(maxen));
      return result;
    }
    if (e0 < e1 && en < e1) {
      result = fcn->Eval(e1) * (en - e0) / (e1 - e0);
      return result;
    }
    result = fcn->Eval(en) / (1. + radFCN(en));
    return result;
  };
  /**
   * Converting the detection efficiency to a form of std::function
   */
  std::function<double(double, double)> eff =
      [teff](double x, double en) {
        int bin = teff->FindFixBin(x, en);
        double result = teff->GetEfficiency(bin);
        return result;
      };
  std::function<double(double)> vcs_fcn_no_spread =
      [opts, &born_fcn, &eff, &teff](double en) {
        if (en <= opts.thsd) {
          return 0.;
        }
        const double s_threshold = opts.thsd * opts.thsd;
        const double s = en * en;
        double result = 0;
        if (teff) {
          result = convolutionKuraevFadin(en, born_fcn, 0, 1. - s_threshold / s, eff);
        } else {
            result = convolutionKuraevFadin(en, born_fcn, 0, 1. - s_threshold / s);
        }
        return result;
      };
  std::function<double(double)> vcs_fcn =
      [vmap, opts, &vcs_fcn_no_spread](double en) {
        double result = 0;
        if (vmap.count("enable-energy-spread")) {
          result = gaussian_conv(en, opts.energy_spread * opts.energy_spread,
                                 vcs_fcn_no_spread);
        } else {
          result = vcs_fcn_no_spread(en);
        }
        return result;
      };
  Eigen::VectorXd tmpRad = Eigen::VectorXd::Zero(opts.n);
  for (std::size_t iter = 0; iter < opts.niter; ++iter) {
    std::cout << "ITER: " << iter << " / " << opts.niter << std::endl;
    for (std::size_t i = 0; i < opts.n; ++i) {
      tmpRad(i) = vcs_fcn(ecm(i)) / born_fcn(ecm(i)) - 1;
      if (std::isnan(tmpRad(i))) {
        tmpRad(i) = 0;
      }
    }
    std::copy(tmpRad.data(), tmpRad.data() + opts.n, radCorr.data());
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
    acc = gsl_interp_accel_alloc ();
    spline = gsl_spline_alloc(gsl_interp_linear, opts.n);
    gsl_spline_init (spline, ecm.data(), radCorr.data(), opts.n);
  }
  std::function<double(double*, double*)> rad_fcn =
      [&radFCN](double* x, double*) {
       double result = radFCN(x[0]);
       return result;
      };
  TF1 rootRadCorrFCN("radCorrFCN", &rad_fcn, minen, maxen, 0);
  rootRadCorrFCN.SetNpx(100);
  std::vector<double> radcorrs;
  std::vector<double> cs;
  std::vector<double> csErr;
  radcorrs.reserve(vcs->GetN());
  cs.reserve(vcs->GetN());
  csErr.reserve(vcs->GetN());
  for (int i = 0; i < vcs->GetN(); ++i) {
    double energy = (vcs->GetX())[i];
    radcorrs.push_back(radFCN(energy));
    cs.push_back((vcs->GetY())[i] /  (1. + radcorrs[i]));
    csErr.push_back((vcs->GetEY())[i] /  (1. + radcorrs[i]));
  }
  TGraph radGraph(vcs->GetN(), vcs->GetX(), radcorrs.data());
  auto bornCS = new TGraphErrors(vcs->GetN(), vcs->GetX(), cs.data(), 0, csErr.data());
  auto ofl = TFile::Open(opts.ofname.c_str(), "recreate");
  ofl->cd();
  rootRadCorrFCN.Write();
  radGraph.Write("radcorr");
  vcs->Write("vcs");
  bornCS->Write("bcs");
  ofl->Close();
  delete fl;
  delete bornCS;
  delete vcs;
  delete fcn;
  if (teff) {
    delete teff;
  }
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  return 0;
}
