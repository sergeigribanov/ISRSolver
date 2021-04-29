#include <boost/program_options.hpp>
#include <iostream>
#include <string>
#include "ISRSolverTikhonov.hpp"
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
   * Regularization parameter
   */
  double lambda;
  /**
   * Name of the visible cross section
   * graph (TGraphErrors)
   */
  std::string vcs_name;
  /**
   * Name of the detection efficiency
   * object (TEfficiency)
   */
  std::string efficiency_name;
  /**
   * Path to the input .root file that contains
   * the visible cross graph (TGraphErrors) and
   * detection efficiency (TEfficiency)
   */
  std::string ifname;
  /**
   * Output file path
   */
  std::string ofname;
  /**
   * Path to the .json file with interpolation
   */
  std::string interp;
} CmdOptions;

/**
 * Setting up program options
 */
void setOptions(po::options_description* desc, CmdOptions* opts) {
  desc->add_options()
      ("help,h", "help message")
      ("thsd,t", po::value<double>(&(opts->thsd)), "threshold energy (GeV)")
      ("enable-energy-spread,g", "enable energy spread")
      ("use-solution-norm2,s",
       "use the following regularizator: lambda*||solution||^2 if this option is enabled, use lambda*||d(solution) / dE||^2 otherwise")
      ("lambda,l", po::value<double>(&(opts->lambda)), "regularization parameter (lambda)")
      ("vcs-name,v", po::value<std::string>(&(opts->vcs_name))->default_value("vcs"),
       "name of the visible cross section graph (TGraphErrors*)")
      ("efficiency-name,e", po::value<std::string>(&(opts->efficiency_name)),
       "name of a detection efficiency object (TEfficiency*)")
      ("ifname,i",  po::value<std::string>(&(opts->ifname))->default_value("vcs.root"),
       "path to input file")
      ("ofname,o", po::value<std::string>(&(opts->ofname))->default_value("bcs.root"),
       "path to output file")
      ("interp,r", po::value<std::string>(&(opts->interp)),
       "path to JSON file with interpolation settings");
}

/**
 * Help message
 */
void help(const po::options_description& desc) {
  std::cout << desc << std::endl;
}

int main(int argc, char* argv[]) {
  po::options_description desc("   Solver using the Tikhonov regularization method. Allowed options:");
  CmdOptions opts;
  setOptions(&desc, &opts);
  po::variables_map vmap;
  po::store(po::parse_command_line(argc, argv, desc), vmap);
  po::notify(vmap);
  if (vmap.count("help")) {
    help(desc);
    return 0;
  }
  /**
   * Creating Tikhonov solver
   */
  ISRSolverTikhonov solver(opts.ifname, {
	    .efficiencyName = opts.efficiency_name,
	    .visibleCSGraphName = opts.vcs_name,
	    .thresholdEnergy = opts.thsd});
  if (vmap.count("enable-energy-spread")) {
    solver.enableEnergySpread();
  }
  if (vmap.count("interp")) {
    solver.setRangeInterpSettings(opts.interp);
  }
  if (vmap.count("lambda")) {
    /**
     * Setting regularization parameter
     */
    solver.setLambda(opts.lambda);
  }
  if (vmap.count("use-solution-norm2")) {
    solver.disableDerivNorm2Regularizator();
  }
  /**
   * Finding a numerical solution
   */
  solver.solve();
  /**
   * Saving result into the output file
   */
  solver.save(opts.ofname,
              {.visibleCSGraphName = opts.vcs_name,
               .bornCSGraphName = "bcs"});
  solver.printConditionNumber();
  return 0;
}
