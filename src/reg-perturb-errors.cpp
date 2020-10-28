#include <boost/program_options.hpp>
#include <iostream>
#include <string>
#include <nlopt.hpp>
#include <TGraph.h>
#include <TGraphErrors.h>
#include "ISRSolverTikhonov.hpp"
#include "utils.hpp"
namespace po = boost::program_options;

typedef struct {
  double thsd;
  double alpha_min;
  double alpha_max;
  int n;
  std::string vcs_name;
  std::string efficiency_name;
  std::string ifname;
  std::string ofname;
  std::string interp;
} CmdOptions;

void setOptions(po::options_description* desc, CmdOptions* opts) {
  desc->add_options()("help,h",
                      "A simple tool designed in order to find numerical"
                      "solution of the Kuraev-Fadin equation.")
    ("thsd,t", po::value<double>(&(opts->thsd)), "Threshold (GeV).")
    ("n-points,n", po::value<int>(&(opts->n))->default_value(100),
     "Number of points in regularization error graph.")
    ("alpha-min,m", po::value<double>(&(opts->alpha_min))->default_value(0.),
     "Minimum value of Tikhonov's regularization parameter.")
    ("alpha-max,x", po::value<double>(&(opts->alpha_max))->default_value(1.),
     "Maximum value of Tikhonov's regularization parameter.")
    ("enable-solution-positivity,p", "Setting positive limits to solution")
    ("disable-solution-norm,f", 
     "Disable solution norm in Tihonov's regularization functional")
    ("disable-solution-derivative-norm,d",
     "Disable solution derivative norm in Tikhonov's regularization functional")
    ("vcs-name,v", po::value<std::string>(&(opts->vcs_name))->default_value("vcs"),
      "Name of the visible cross section graph.")
    ("efficiency-name,e", po::value<std::string>(&(opts->efficiency_name)),
     "TEfficiency object name")
    ( "ifname,i",
      po::value<std::string>(&(opts->ifname))->default_value("vcs.root"),
      "Path to input file.")(
      "ofname,o",
      po::value<std::string>(&(opts->ofname))->default_value("bcs.root"),
      "Path to output file.")("interp,r",
                              po::value<std::string>(&(opts->interp)),
                              "Path to JSON file with interpolation settings.");
}

void help(const po::options_description& desc) {
  std::cout << desc << std::endl;
}

double alphaObjective(unsigned n, const double* palpha, double* grad, void* solver) {
   auto sp = reinterpret_cast<ISRSolverTikhonov*>(solver);
   sp->setAlpha(*palpha);
   sp->solve();
   return sp->evalCurvature();
}

int main(int argc, char* argv[]) {
  po::options_description desc("Allowed options:");
  CmdOptions opts;
  setOptions(&desc, &opts);
  po::variables_map vmap;
  po::store(po::parse_command_line(argc, argv, desc), vmap);
  po::notify(vmap);
  if (vmap.count("help")) {
    help(desc);
    return 0;
  }

  TGraphErrors* bcsSRC;
  TGraphErrors* vcsSRC;
  TGraphErrors* vcsP;
  auto ifl = TFile::Open(opts.ifname.c_str(), "read");
  bcsSRC = dynamic_cast<TGraphErrors*>(ifl->Get("bcsSRC")->Clone("src_bcs"));
  vcsSRC = dynamic_cast<TGraphErrors*>(ifl->Get("vcsSRC")->Clone("src_vcs"));
  vcsP = dynamic_cast<TGraphErrors*>(ifl->Get("vcs")->Clone("src_perturbated"));
  ifl->Close();
  delete ifl;
  Eigen::VectorXd bcsOrig = Eigen::VectorXd::Zero(bcsSRC->GetN());
  Eigen::VectorXd vcsPerturbation = Eigen::VectorXd::Zero(bcsSRC->GetN());
  for (int i = 0; i < bcsSRC->GetN(); ++i) {
    bcsOrig(i) = bcsSRC->GetY()[i];
    vcsPerturbation(i) = vcsP->GetY()[i] - vcsSRC->GetY()[i];
  }

  delete bcsSRC;
  delete vcsSRC;
  delete vcsP;
  
  auto solver = new ISRSolverTikhonov(opts.ifname, {
      .efficiencyName = opts.efficiency_name,
      .visibleCSGraphName = opts.vcs_name,
      .thresholdEnergy = opts.thsd,
      .energyUnitMeVs = false});
  
  if (vmap.count("interp")) {
    solver->setInterpSettings(opts.interp);
  }
  if (vmap.count("disable-solution-derivative-norm")) {
    solver->disableSolutionDerivativeNorm2();
  }
  if (vmap.count("enable-solution-positivity")) {
    solver->enableSolutionPositivity();
  }
  TGraph relRegErr;
  TGraph relPerturpErr;
  double alphaStep = (opts.alpha_max - opts.alpha_min) / (opts.n - 1);
  for (int i = 0; i < opts.n; ++i) {
    double alpha = opts.alpha_min + i * alphaStep;
    solver->setAlpha(alpha);
    solver->solve();
    relRegErr.SetPoint(i, alpha, solver->evalApproxRegRelativeError(bcsOrig));
    relPerturpErr.SetPoint(i, alpha, solver->evalApproxPerturbRelativeError(bcsOrig, vcsPerturbation));
    std::cout << "progress: " << i + 1 << " / " << opts.n << std::endl;
    std::cout << "alpha = " << alpha << std::endl;
  }
  delete solver;
  auto fl = TFile::Open(opts.ofname.c_str(), "recreate");
  fl->cd();
  relRegErr.Write("relRegErr");
  relPerturpErr.Write("relPerturpErr");
  fl->Close();
  delete fl;
  return 0;
}
