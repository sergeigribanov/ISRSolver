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
      ("help,h", "help message")
      ("thsd,t", po::value<double>(&(opts->thsd)), "threshold (GeV)")
      ("enable-energy-spread,g", "enable energy spread")
      ("ampl,l", po::value<double>(&(opts->ampl))->default_value(1.e+4),
       "initial chi-square amplitude")
      ("num-rnd-draws,n", po::value<int>(&(opts->n))->default_value(100),
       "number of visible cross section random draws")
      ("vcs-name,v",
       po::value<std::string>(&(opts->vcs_name))->default_value("vcs"),
       "name of the visible cross section graph (TGraphErrors*)")
      ("efficiency-name,e", po::value<std::string>(&(opts->efficiency_name)),
       "name of the detection efficiency object (TEfficiency*)")
      ("use-model,u", po::value<std::string>(&(opts->path_to_model)),
       "path to the file with the model Born and visible cross section in the form of graphs (TGraphErrors*)")
      ("model-bcs-name,b",
       po::value<std::string>(&(opts->name_of_model_bcs))->default_value("bcsSRC"),
       "name of the model Born cross section graph (TGraphErrors*)")
      ("model-vcs-name,c",
       po::value<std::string>(&(opts->name_of_model_vcs))->default_value("vcsBlured"),
       "name of the model visible cross graph (TGraphErrors*)")
      ( "ifname,i",
        po::value<std::string>(&(opts->ifname))->default_value("vcs.root"),
        "path to input file")
      ("ofname,o",
       po::value<std::string>(&(opts->ofname))->default_value("isrsolver-chi2-distribution.root"),
       "path to output file")
      ("interp,r",
       po::value<std::string>(&(opts->interp)),
       "path to JSON file with interpolation settings");
}

void help(const po::options_description& desc) {
  std::cout << desc << std::endl;
}

int main(int argc, char* argv[]) {
  po::options_description desc("   This tool calculates chi-square in multiple numerical experiments using "
                               "the naive method and fills chi-square histogram. Allowed options");
  CmdOptions opts;
  setOptions(&desc, &opts);
  po::variables_map vmap;
  po::store(po::parse_command_line(argc, argv, desc), vmap);
  po::notify(vmap);
  if (vmap.count("help")) {
    help(desc);
    return 0;
  }
  ISRSolverSLE solver(
      opts.ifname,
      {.efficiencyName = opts.efficiency_name,
       .visibleCSGraphName = opts.vcs_name,
       .thresholdEnergy = opts.thsd});
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
