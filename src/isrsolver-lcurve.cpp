#include <boost/program_options.hpp>
#include <iostream>
#include <string>
#include <nlopt.hpp>
#include "ISRSolverTikhonov.hpp"
#include "utils.hpp"
namespace po = boost::program_options;

typedef struct {
  double thsd;
  double alpha;
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
      ("enable-energy-spread,g", "Enable energy spread")
      // ("enable-solution-positivity,p", "Setting positive limits to solution")
      ("disable-solution-norm,f", 
       "Disable solution norm in Tihonov's regularization functional")
      ("disable-solution-derivative-norm,d",
       "Disable solution derivative norm in Tikhonov's regularization functional")
      ("alpha,a", po::value<double>(&(opts->alpha))->default_value(1.e-9),
       "Thikhonov's regularization parameter.")(
           "vcs-name,v", po::value<std::string>(&(opts->vcs_name))->default_value("vcs"),
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
   if (grad) {
     grad[0] = sp->evalLCurveCurvatureDerivative();
   }
   return sp->evalLCurveCurvature();
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
  auto solver = new ISRSolverTikhonov(opts.ifname, {
      .efficiencyName = opts.efficiency_name,
      .visibleCSGraphName = opts.vcs_name,
      .thresholdEnergy = opts.thsd,
      .energyUnitMeVs = false});
  std::vector<double> z(1, opts.alpha);
  if (vmap.count("interp")) {
    solver->setInterpSettings(opts.interp);
  }
  if (vmap.count("disable-solution-derivative-norm")) {
    solver->disableSolutionDerivativeNorm2();
  }
  if (vmap.count("enable-energy-spread")) {
    solver->enableEnergySpread();
  }
  nlopt::opt opt(nlopt::LD_MMA, 1);
  std::vector<double> lowerBounds(1, 0);
  opt.set_lower_bounds(lowerBounds);
  opt.set_min_objective(alphaObjective, solver);
  opt.set_xtol_rel(1.e-6);
  double minf;
  opt.optimize(z, minf);
  std::cout << "alpha = " << z[0] << std::endl;
  std::cout << "curvature = " << minf << std::endl;
  solver->save(opts.ofname,
               {.visibleCSGraphName = opts.vcs_name, .bornCSGraphName = "bcs"});
  delete solver;
  return 0;
}
