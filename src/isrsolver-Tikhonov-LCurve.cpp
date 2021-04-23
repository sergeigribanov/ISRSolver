#include <iostream>
#include <string>
#include <boost/program_options.hpp>
#include <nlopt.hpp>
#include "ISRSolverTikhonov.hpp"
namespace po = boost::program_options;

typedef struct {
  double thsd;
  double lambda;
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
      ("use-solution-norm2,s",
       "regularization: lambda*||solution||^2 if enabled, lambda*||d(solution) / dE||^2 otherwise.")
      ("lambda,l", po::value<double>(&(opts->lambda))->default_value(1.e-9),
       "regularization parameter")(
           "vcs-name,v", po::value<std::string>(&(opts->vcs_name))->default_value("vcs"),
           "name of the visible cross section graph (TGraphErrors*)")
      ("efficiency-name,e", po::value<std::string>(&(opts->efficiency_name)),
       "name of the detection efficiency object (TEfficiency*)")
      ( "ifname,i",
        po::value<std::string>(&(opts->ifname))->default_value("vcs.root"),
        "path to input file")
      ("ofname,o",
       po::value<std::string>(&(opts->ofname))->default_value("bcs.root"),
       "path to output file")
      ("interp,r",
       po::value<std::string>(&(opts->interp)),
       "path to JSON file with interpolation settings");
}

void help(const po::options_description& desc) {
  std::cout << desc << std::endl;
}

double lambdaObjective(unsigned n, const double* plambda, double* grad, void* solver) {
   auto sp = reinterpret_cast<ISRSolverTikhonov*>(solver);
   sp->setLambda(*plambda);
   sp->solve();
   if (grad) {
     grad[0] = sp->evalLCurveCurvatureDerivative();
   }
   return sp->evalLCurveCurvature();
}

int main(int argc, char* argv[]) {
  po::options_description desc("   Solver that uses Tikhonov regularization in combination with an L-curve criterion. Allowed options:");
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
      .thresholdEnergy = opts.thsd});
  std::vector<double> z(1, opts.lambda);
  if (vmap.count("interp")) {
    solver->setRangeInterpSettings(opts.interp);
  }
  if (vmap.count("use-solution-norm2")) {
    solver->disableDerivNorm2Regularizator();
  }
  if (vmap.count("enable-energy-spread")) {
    solver->enableEnergySpread();
  }
  nlopt::opt opt(nlopt::LD_MMA, 1);
  std::vector<double> lowerBounds(1, 0);
  opt.set_lower_bounds(lowerBounds);
  opt.set_min_objective(lambdaObjective, solver);
  opt.set_xtol_rel(1.e-6);
  double minf;
  opt.optimize(z, minf);
  std::cout << "lambda = " << z[0] << std::endl;
  std::cout << "curvature = " << minf << std::endl;
  solver->save(opts.ofname,
               {.visibleCSGraphName = opts.vcs_name, .bornCSGraphName = "bcs"});
  delete solver;
  return 0;
}
