#include <iostream>
#include <string>
#include <boost/program_options.hpp>
#include "IterISRInterpSolver.hpp"
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
   * Threshold energy
   */
  double thsd;
  /**
   * Name of the visible cross section graph (TGraphErrors)
   */
  std::string vcs_name;
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
  /**
   * Path to the file with interpolation settings
   */
  std::string interp;
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
      ("thsd,t", po::value<double>(&(opts->thsd)), "threshold (GeV)")
      ("vcs-name,v", po::value<std::string>(&(opts->vcs_name))->default_value("vcs"),
       "name of the visible cross section graph (TGraphErrors*)")
      ("efficiency-name,e", po::value<std::string>(&(opts->efficiency_name)),
       "name of the detection efficiency object (TEfficiency*)")
      ( "ifname,i",
        po::value<std::string>(&(opts->ifname))->default_value("vcs.root"),
        "Path to input file.")
      ("ofname,o",
       po::value<std::string>(&(opts->ofname))->default_value("bcs.root"),
       "Path to output file.")
      ("interp,r", po::value<std::string>(&(opts->interp)),
       "Path to JSON file with interpolation settings.");
}

/**
 * Help message
 */
void help(const po::options_description& desc) {
  std::cout << desc << std::endl;
}

int main(int argc, char* argv[]) {
  po::options_description desc("   Common used iterative solver. Allowed options");
  CmdOptions opts;
  setOptions(&desc, &opts);
  po::variables_map vmap;
  po::store(po::parse_command_line(argc, argv, desc), vmap);
  po::notify(vmap);
  if (vmap.count("help")) {
    help(desc);
    return 0;
  }
  IterISRInterpSolver solver(opts.ifname, {
      .efficiencyName = opts.efficiency_name,
      .visibleCSGraphName = opts.vcs_name,
      .thresholdEnergy = opts.thsd});
  solver.setNumOfIters(opts.niter);
  if (vmap.count("enable-energy-spread")) {
    solver.enableEnergySpread();
  }
  solver.solve();
  solver.save(opts.ofname,
               {.visibleCSGraphName = opts.vcs_name, .bornCSGraphName = "bcs"});
  return 0;
}
