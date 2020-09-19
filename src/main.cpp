#include <TAxis.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>

#include <boost/program_options.hpp>
#include <cstdlib>
#include <iostream>
#include <string>

#include "ISRSolver.hpp"
#include "utils.hpp"
namespace po = boost::program_options;

typedef struct {
  double thsd;
  double spoint;
  std::string gname;
  std::string lbcs;
  std::string ifname;
  std::string ofname;
} CmdOptions;

void setOptions(po::options_description* desc, CmdOptions* opts) {
  desc->add_options()("help",
                      "A simple tool designed in order to find numerical"
                      "solution of the Kuraev-Fadin equation.")(
      "thsd", po::value<double>(&(opts->thsd)), "Threshold (GeV).")(
      "gname", po::value<std::string>(&(opts->gname))->default_value("mcs"),
      "Name of the measured cross section graph.")(
      "ifname",
      po::value<std::string>(&(opts->ifname))->default_value("mcs.root"),
      "Path to input file.")(
      "ofname",
      po::value<std::string>(&(opts->ofname))->default_value("bcs.root"),
      "Path to output file.");
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
  ISRSolver solver(opts.ifname, {.measuredCSGraphName = opts.gname,
                                 .thresholdEnergy = opts.thsd,
                                 .energyUnitMeVs = false});
  // std::vector<InterpPtSettings> interpSettings(30, {.numCoeffs = 5,
  // .nPointsLeft = 2}); interpSettings[0] = {.numCoeffs = 2, .nPointsLeft = 1};
  // interpSettings[28] = {.numCoeffs = 3, .nPointsLeft = 1};
  // interpSettings[29] = {.numCoeffs = 2, .nPointsLeft = 1};
  // solver.setInterpSettings(interpSettings);
  solver.solve();
  solver.save(opts.ofname,
              {.measuredCSGraphName = opts.gname, .bornCSGraphName = "bcs"});
  return 0;
}
