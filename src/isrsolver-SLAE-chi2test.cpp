#include <iostream>
#include <boost/program_options.hpp>
#include "Chi2Test.hpp"
namespace po = boost::program_options;

typedef struct {
  double thsd;
  int n;
  double ampl;
  std::string path_to_model;
  std::string name_of_model_bcs;
  std::string name_of_model_vcs;
  std::string vcs_name;
  std::string efficiency_name;
  std::string ifname;
  std::string ofname;
  std::string interp;
} CmdOptions;

void setOptions(po::options_description* desc, CmdOptions* opts) {
  desc->add_options()
      ("help,h", "Help.")
      ("thsd,t", po::value<double>(&(opts->thsd)), "Threshold (GeV).")
      ("enable-energy-spread,g", "Enable energy spread")
      ("ampl,l", po::value<double>(&(opts->ampl))->default_value(1.e+4),
       "Initial chi-square amplitude.")
      ("num-rnd-draws,n", po::value<int>(&(opts->n))->default_value(100),
       "Number of visible cross section random draws.")
      ("vcs-name,v",
       po::value<std::string>(&(opts->vcs_name))->default_value("vcs"),
       "Name of the visible cross section graph.")
      ("efficiency-name,e", po::value<std::string>(&(opts->efficiency_name)),
       "TEfficiency object name.")
      ("use-model,u", po::value<std::string>(&(opts->path_to_model)),
       "Path to the file with the model Born and visible cross section TGraphErrors (if needed).")
      ("model-bcs-name,b",
       po::value<std::string>(&(opts->name_of_model_bcs))->default_value("bcsSRC"),
       "Name of the model Born cross section TGraphErrors function")
      ("model-vcs-name,c",
       po::value<std::string>(&(opts->name_of_model_vcs))->default_value("vcsBlured"),
       "Name of the model visible cross section TGraphErrors function")
      ( "ifname,i",
        po::value<std::string>(&(opts->ifname))->default_value("vcs.root"),
        "Path to input file.")
      ("ofname,o",
       po::value<std::string>(&(opts->ofname))->default_value("isrsolver-chi2-distribution.root"),
       "Path to output file.")
      ("interp,r",
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
  ISRSolverSLAE solver(
      opts.ifname,
      {.efficiencyName = opts.efficiency_name,
       .visibleCSGraphName = opts.vcs_name,
       .thresholdEnergy = opts.thsd,
       .energyUnitMeVs = false});
  if (vmap.count("interp")) {
    solver.setRangeInterpSettings(opts.interp);
  }
  if (vmap.count("enable-energy-spread")) {
    solver.enableEnergySpread();
  }
  if (vmap.count("use-model")) {
    chi2TestModel(opts.n, opts.ampl, &solver,
                  opts.path_to_model,
                  opts.name_of_model_vcs,
                  opts.name_of_model_bcs,
                  opts.ofname);
  } else {
    chi2TestData(opts.n, opts.ampl, &solver, opts.ofname);
  }
  return 0;
}
