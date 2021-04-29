#include <string>
#include <boost/program_options.hpp>
#include "RatioTest.hpp"
namespace po = boost::program_options;

/**
 * A part of program options
 */
typedef struct {
  /**
   * Threshold energy
   */
  double thsd;
  /**
   * Number of points in the ratio graph
   */
  int n;
  /**
   * Path to the .root file that contains
   * the model Born cross section graph
   * (TGraphErrors)
   */
  std::string path_to_model;
  /**
   * Name of the model Born cross section
   * graph (TGraphErrors)
   */
  std::string name_of_model_bcs;
  /**
   * Name of the visible Born cross section
   * graph (TGraphErrors)
   */
  std::string name_of_model_vcs;
  /**
   * Name of the visible cross section
   * graph (TGraphErrors)
   */
  std::string vcs_name;
  /**
   * Name of the detection efficiency object
   * (TEfficiency)
   */
  std::string efficiency_name;
  /**
   * Path to the input .root file
   * that contains the visible cross
   * section graph and the detection
   * efficiency
   */
  std::string ifname;
  /**
   * Path to the output file
   */
  std::string ofname;
  /**
   * Path to the .json file with
   * interpolation settings
   */
  std::string interp;
} CmdOptions;

/**
 * Setting up program options
 */
void setOptions(po::options_description* desc, CmdOptions* opts) {
  desc->add_options()
      ("help,h", "help message")
      ("thsd,t", po::value<double>(&(opts->thsd)), "threshold (GeV)")
      ("num-rnd-draws,n", po::value<int>(&(opts->n))->default_value(100),
       "number of visible cross section random draws")
      ("enable-energy-spread,g", "enable energy spread")
      ("vcs-name,v",
       po::value<std::string>(&(opts->vcs_name))->default_value("vcs"),
       "name of the visible cross section graph (TGraphErrors*)")
      ("efficiency-name,e", po::value<std::string>(&(opts->efficiency_name)),
       "name of the detection efficiency object (TEfficiency*)")
      ("use-model,u", po::value<std::string>(&(opts->path_to_model)),
       "path to the file with the model Born and visible cross sections in the form of graphs (TGraphErrors*)")
      ("model-bcs-name,b",
       po::value<std::string>(&(opts->name_of_model_bcs))->default_value("bcsSRC"),
       "name of the model Born cross section graph (TGraphErrors*)")
      ("model-vcs-name,c",
       po::value<std::string>(&(opts->name_of_model_vcs))->default_value("vcsBlured"),
       "name of the model visible cross section graph (TGraphErrors*)")
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

/**
 * Help message
 */
void help(const po::options_description& desc) {
  std::cout << desc << std::endl;
}

int main(int argc, char* argv[]) {
  po::options_description desc(
      "   This tool is calculates ratio between a numerical solution and a model Born cross section. "
      "This ratio is averaged over multiple numerical experiments. The naive method is used in each "
      "numerical experiment to obtain the numerical solution for the Born cross section. Allowed options");
  CmdOptions opts;
  setOptions(&desc, &opts);
  po::variables_map vmap;
  po::store(po::parse_command_line(argc, argv, desc), vmap);
  po::notify(vmap);
  if (vmap.count("help")) {
    help(desc);
    return 0;
  }
  if (!vmap.count("use-model")) {
    std::cout << "[!] use-model command option is obligatory" << std::endl;
    return 0;
  }
  /**
   * Creating SLE solver
   */
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
  /**
   * Obtaining mean ratio of the numerical solution
   * to the model Born cross section
   */
  ratioTestModel(&solver,
                 {.n = opts.n,
                  .modelPath = opts.path_to_model,
                  .modelVCSName = opts.name_of_model_vcs,
                  .modelBCSName = opts.name_of_model_bcs,
                  .outputPath = opts.ofname});
  return 0;
}
