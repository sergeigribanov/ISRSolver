#include "ISRSolverVCSFitter.hpp"
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnMinos.h>
#include <Minuit2/MnUserParameters.h>
#include <Minuit2/MnUserParameterState.h>
#include <Minuit2/FunctionMinimum.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <boost/program_options.hpp>
#include <iostream>
#include <vector>
namespace po = boost::program_options;

using namespace ROOT::Minuit2;

typedef struct {
  double thsd;
  std::string vcs_name;
  std::string ifname;
  std::string ofname;
} CmdOptions;

void setOptions(po::options_description *desc, CmdOptions *opts) {
  desc->add_options()("help,h", "help message")(
      "thsd,t", po::value<double>(&(opts->thsd)), "threshold energy (GeV)")(
      "enable-energy-spread,g", "enable energy spread")(
      "vcs-name,v",
      po::value<std::string>(&(opts->vcs_name))->default_value("vcs"),
      "name of a visible cross section graph (TGraphErrors*)")(
      "ifname,i",
      po::value<std::string>(&(opts->ifname))->default_value("vcs.root"),
      "path to input file")(
      "ofname,o",
      po::value<std::string>(&(opts->ofname))->default_value("bcs.root"),
      "path to output file");
}

void help(const po::options_description &desc) {
  std::cout << desc << std::endl;
}

Double_t fitstepbw(double x, const std::vector<double>& par) {
  // for (auto el : par)
  //   std::cout << el << " ";
  // std::cout << std::endl;

  x *= 1.e+3;
  Double_t arg1 = 0;
  Double_t bin = 5;
  Double_t thr = 1390.;
  Double_t r = x; //-thr;

  Double_t M = par[4];
  Double_t W = par[5];
  Double_t Vmax = par[3]; //*bin*2./W/3.1415927;
  Double_t BW = (W * W / 4.) * Vmax / ((x - M) * (x - M) + W * W / 4.);

  Double_t M1 = par[7];
  Double_t W1 = par[8];
  Double_t V1max = par[6]; //*bin*2./W/3.1415927;
  Double_t BW1 =
      (W1 * W1 / 4.) * V1max / ((x - M1) * (x - M1) + W1 * W1 / 4.);

  if (par[2] != 0)
    arg1 = (x - par[1]) / (par[2]);

  Double_t FS = ((1 - x * x / thr / thr) / (1 - M * M / thr / thr));

  Double_t fitval;

  fitval = FS * (par[9] + par[10] * r * r * r * r +
                 par[0] / (1. + TMath::Exp(arg1)) + BW + BW1);

  if (r < thr)
    fitval = 0;
  if (fitval < 0)
    fitval = 0;

  // std::cout << "x = " << x << ", fitval = " << fitval << std::endl;
  // std::cout << "_________" << std::endl;

  return fitval;
}

int main(int argc, char *argv[]) {
  po::options_description desc("  Visible cross-section fitter");
  CmdOptions opts;
  setOptions(&desc, &opts);
  po::variables_map vmap;
  po::store(po::parse_command_line(argc, argv, desc), vmap);
  po::notify(vmap);
  if (vmap.count("help")) {
    help(desc);
    return 0;
  }
  auto ifl = TFile::Open(opts.ifname.c_str(), "read");
  auto vcs = dynamic_cast<TGraphErrors*>(ifl->Get(opts.vcs_name.c_str()));
  std::size_t n = vcs->GetN();
  std::vector<double> x(n);
  std::vector<double> y(n);
  std::vector<double> ex(n);
  std::vector<double> ey(n);
  std::copy(vcs->GetX(), vcs->GetX() + n, x.begin());
  std::copy(vcs->GetY(), vcs->GetY() + n, y.begin());
  std::copy(vcs->GetEX(), vcs->GetEX() + n, ex.begin());
  std::copy(vcs->GetEY(), vcs->GetEY() + n, ey.begin());

  ifl->Close();
  delete ifl;
  ISRSolverVCSFitFunction fitter(n, opts.thsd, x.data(),
                                 y.data(), ex.data(),
                                 ey.data(), fitstepbw);
  if (vmap.count("enable-energy-spread")) {
    fitter.enableEnergySpread();
  }
  MnUserParameters upar;
  upar.Add("0", 0.25, 1.e-5);
  upar.Add("1", 1877., 1.e-5); // , 1875., 1880.
  upar.Add("2", 2.9, 1.e-5);   // , 0.1, 3.
  upar.Add("3", 0.1, 1.e-5);   // , 0., 1.0
  upar.Add("4", 1800., 1.e-5); // , 1750., 2000.
  upar.Add("5", 200., 1.e-5);  // ,  100., 300.
  upar.Add("6", 0.1, 1.e-5);   // , 0., 1.
  upar.Add("7", 1720., 1.e-5); // , 1650., 1750.
  upar.Add("8", 150., 1.e-5);  // 5, 50., 300.
  upar.Add("9", 0., 1.e-5);
  upar.Add("10", 0., 1.e-5);

  MnMigrad migrad(fitter, upar);
  FunctionMinimum min = migrad();
  std::cout << "minimum: " << min << std::endl;
  fitter.saveResults(min, opts.ofname);
  return 0;
}
