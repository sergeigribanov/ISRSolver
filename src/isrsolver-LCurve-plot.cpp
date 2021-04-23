#include <cmath>
#include <iostream>
#include <string>
#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <TGraph.h>
#include <TFile.h>
#include "ISRSolverTikhonov.hpp"
namespace po = boost::program_options;

typedef struct {
  double thsd;
  double lambda_min;
  double lambda_max;
  int lambda_n;
  std::string vcs_name;
  std::string efficiency_name;
  std::string ifname;
  std::string ofname;
  std::string interp;
  std::string solver;
} CmdOptions;

void setOptions(po::options_description* desc, CmdOptions* opts) {
  desc->add_options()
      ("help,h", "help message")
      ("thsd,t", po::value<double>(&(opts->thsd)), "threshold energy (GeV)")
      ("enable-energy-spread,g", "enable energy spread")
      ("use-solution-norm2,s",
       "use the following regularizator: lambda*||solution||^2 if enabled, use lambda*||d(solution) / dE||^2 otherwise")
      ("lambda-min,m", po::value<double>(&(opts->lambda_min))->default_value(1.e-9),
       "minimum value of regularization parameter")
      ("lambda-max,x", po::value<double>(&(opts->lambda_max))->default_value(1.0),
       "maximum value of regularization parameter.")
      ("lambda-n,n", po::value<int>(&(opts->lambda_n))->default_value(10),
       "number of steps in regularization parameter.")
      ("vcs-name,v", po::value<std::string>(&(opts->vcs_name))->default_value("vcs"),
       "name of a visible cross section graph (TGraphErrors*)")
      ("efficiency-name,e", po::value<std::string>(&(opts->efficiency_name)),
       "name of a detection efficiency object (TEfficiency*)")
      ( "ifname,i",
        po::value<std::string>(&(opts->ifname))->default_value("vcs.root"),
        "path to input file")
      ( "ofname,o",
        po::value<std::string>(&(opts->ofname))->default_value("bcs.root"),
        "path to output file")
      ("interp,r",
       po::value<std::string>(&(opts->interp)),
       "path to JSON file with interpolation settings");
}

void help(const po::options_description& desc) {
  std::cout << desc << std::endl;
}

int main(int argc, char* argv[]) {
  po::options_description desc("Drawing L-Curve. Allowed options:");
  CmdOptions opts;
  setOptions(&desc, &opts);
  po::variables_map vmap;
  po::store(po::parse_command_line(argc, argv, desc), vmap);
  po::notify(vmap);
  if (vmap.count("help")) {
    help(desc);
    return 0;
  }
  ISRSolverTikhonov solver(opts.ifname, {
      .efficiencyName = opts.efficiency_name,
      .visibleCSGraphName = opts.vcs_name,
      .thresholdEnergy = opts.thsd});
  if (vmap.count("interp")) {
    solver.setRangeInterpSettings(opts.interp);
  }
  if (vmap.count("use-solution-norm2")) {
    solver.disableDerivNorm2Regularizator();
  }
  if (vmap.count("enable-energy-spread")) {
    solver.enableEnergySpread();
  }
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> a;
  std::vector<double> curv;
  x.reserve(opts.lambda_n);
  y.reserve(opts.lambda_n);
  a.reserve(opts.lambda_n);
  curv.reserve(opts.lambda_n);
  double h = std::pow(opts.lambda_max / opts.lambda_min,  1. / (opts.lambda_n - 1));
  for (int i = 0; i < opts.lambda_n; ++i) {
    double lambda = opts.lambda_min * std::pow(h, i);
    std::cout << "[" << i + 1 << "/" << opts.lambda_n << "]" << std::endl;
    std::cout << "lambda = " << lambda << std::endl;
    std::cout << "--------" << std::endl;
    solver.setLambda(lambda);
    solver.solve();
    x.push_back(solver.evalEqNorm2());
    y.push_back(solver.evalSmoothnessConstraintNorm2());
    a.push_back(lambda);
    curv.push_back(solver.evalLCurveCurvature());
  }
  TGraph chi2(opts.lambda_n, a.data(), x.data());
  TGraph lcurve(opts.lambda_n, x.data(), y.data());
  TGraph curvature(opts.lambda_n, a.data(), curv.data());
  auto fl = TFile::Open(opts.ofname.c_str(), "recreate");
  fl->cd();
  chi2.Write("chi2");
  lcurve.Write("lcurve");
  curvature.Write("curvature");
  fl->Close();
  delete fl;
  return 0;
}
