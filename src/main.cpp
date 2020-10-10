#include <TAxis.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>

#include <boost/program_options.hpp>
#include <cstdlib>
#include <iostream>
#include <string>

#include "ISRSolverSLAE.hpp"
#include "ISRSolverTikhonov.hpp"
#include "utils.hpp"
namespace po = boost::program_options;

typedef struct {
  double thsd;
  double alpha;
  double spoint;
  std::string gname;
  std::string lbcs;
  std::string ifname;
  std::string ofname;
  std::string interp;
  std::string solver;
} CmdOptions;

void setOptions(po::options_description* desc, CmdOptions* opts) {
  desc->add_options()("help,h",
                      "A simple tool designed in order to find numerical"
                      "solution of the Kuraev-Fadin equation.")
    ("thsd,t", po::value<double>(&(opts->thsd)), "Threshold (GeV).")
    ("enable-solution-positivity,p", "Setting positive limits to solution")
    ("disable-solution-norm,f", 
     "Disable solution norm in Tihonov's regularization functional")
    ("disable-solution-derivative-norm,d",
     "Disable solution derivative norm in Tikhonov's regularization functional")
    (
      "alpha,a", po::value<double>(&(opts->alpha)),
      "Thikhonov's regularization parameter.")(
      "solver,s", po::value<std::string>(&(opts->solver)),
      "Solver: SLAE, Tikhonov")(
      "gname,g", po::value<std::string>(&(opts->gname))->default_value("vcs"),
      "Name of the measured cross section graph.")(
      "ifname,i",
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
  BaseISRSolver* solver = nullptr;
  if (opts.solver == "SLAE") {
    solver = new ISRSolverSLAE(opts.ifname, {.measuredCSGraphName = opts.gname,
                                             .thresholdEnergy = opts.thsd,
                                             .energyUnitMeVs = false});
  } else if (opts.solver == "Tikhonov") {
    solver =
        new ISRSolverTikhonov(opts.ifname, {.measuredCSGraphName = opts.gname,
                                            .thresholdEnergy = opts.thsd,
                                            .energyUnitMeVs = false});
  }
  if (!solver) {
    std::cerr << "[!] Solver is not set." << std::endl;
    return 1;
  }
  ISRSolverSLAE* solverSLAE = dynamic_cast<ISRSolverSLAE*>(solver);
  if (vmap.count("interp") && solverSLAE) {
    solverSLAE->setInterpSettings(opts.interp);
  }
  ISRSolverTikhonov* solverTikhonov = dynamic_cast<ISRSolverTikhonov*>(solver);
  if (vmap.count("alpha") && solverTikhonov) {
    solverTikhonov->setAlpha(opts.alpha);
  }
  if (vmap.count("disable-solution-norm") && solverTikhonov) {
    solverTikhonov->disableSolutionNorm2();
  }
  if (vmap.count("disable-solution-derivative-norm") && solverTikhonov) {
    solverTikhonov->disableSolutionDerivativeNorm2();
  }
  if (vmap.count("enable-solution-positivity") && solverTikhonov) {
    solverTikhonov->enableSolutionPositivity();
  }
  solver->solve();
  solver->save(opts.ofname,
               {.measuredCSGraphName = opts.gname, .bornCSGraphName = "bcs"});
  delete solver;
  return 0;
}
