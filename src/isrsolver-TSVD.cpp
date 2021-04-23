#include <boost/program_options.hpp>
#include <iostream>
#include <string>
#include "ISRSolverTSVD.hpp"
namespace po = boost::program_options;

typedef struct {
  double thsd;
  int k;
  std::string vcs_name;
  std::string efficiency_name;
  std::string ifname;
  std::string ofname;
  std::string interp;
} CmdOptions;

void setOptions(po::options_description* desc, CmdOptions* opts) {
  desc->add_options()
      ("help,h", "help message")
      ("thsd,t", po::value<double>(&(opts->thsd)), "threshold energy (GeV)")
      ("enable-energy-spread,g", "enable energy spread")
      ("upper-tsvd-index,k", po::value<int>(&(opts->k))->default_value(1), "upper TSVD index")
      ("keep-one,z", "keep only k-th SVD harmonic")
      ("vcs-name,v", po::value<std::string>(&(opts->vcs_name))->default_value("vcs"),
       "name of a visible cross section graph (TGraphErrors*)")
      ("ifname,i", po::value<std::string>(&(opts->ifname))->default_value("vcs.root"),
       "path to input file")
      ("ofname,o", po::value<std::string>(&(opts->ofname))->default_value("bcs.root"),
       "path to output file")
      ("efficiency-name,e", po::value<std::string>(&(opts->efficiency_name)),
       "name of a detection efficiency object (TEfficiency*)")
      ("interp,r", po::value<std::string>(&(opts->interp)),
       "path to JSON file with interpolation settings");
}

void help(const po::options_description& desc) {
  std::cout << desc << std::endl;
}

int main(int argc, char* argv[]) {
  po::options_description desc("   Solver that uses Truncated Singular Value Decomposition (TSVD). Allowed options");
  CmdOptions opts;
  setOptions(&desc, &opts);
  po::variables_map vmap;
  po::store(po::parse_command_line(argc, argv, desc), vmap);
  po::notify(vmap);
  if (vmap.count("help")) {
    help(desc);
    return 0;
  }
  ISRSolverTSVD solver(
      opts.ifname,
      {.efficiencyName = opts.efficiency_name,
       .visibleCSGraphName = opts.vcs_name,
       .thresholdEnergy = opts.thsd}, 1);
  if (vmap.count("enable-energy-spread")) {
    solver.enableEnergySpread();
  }
  if (vmap.count("interp")) {
    solver.setRangeInterpSettings(opts.interp);
  }
  if (vmap.count("upper-tsvd-index")) {
    solver.setUpperTSVDIndex(opts.k);
  }
  if (vmap.count("keep-one")) {
    solver.enableKeepOne();
  }
  solver.solve();
  solver.save(opts.ofname,
              {.visibleCSGraphName = opts.vcs_name,
               .bornCSGraphName = "bcs"});
  solver.printConditionNumber();
}
