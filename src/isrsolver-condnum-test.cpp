#include <boost/program_options.hpp>
#include <iostream>
#include <string>
#include <Eigen/Core>
#include <Eigen/SVD>
#include <TGraph.h>
#include <TFile.h>
#include "ISRSolverSLAE.hpp"
#include "utils.hpp"
namespace po = boost::program_options;

typedef struct {
  double thsd;
  int n;
  double ensigma_max;
  std::string vcs_name;
  std::string efficiency_name;
  std::string ifname;
  std::string ofname;
  std::string interp;
} CmdOptions;

void setOptions(po::options_description* desc, CmdOptions* opts) {
  desc->add_options()
      ("help,h", "help message")
      ("num,n", po::value<int>(&(opts->n))->default_value(10), "number of points")
       ("ensigma-max,m",
        po::value<double>(&(opts->ensigma_max))->default_value(0.1),
        "max center-of-mass energy sigma")
       ("thsd,t", po::value<double>(&(opts->thsd)), "threshold (GeV)")
      ("vcs-name,v", po::value<std::string>(&(opts->vcs_name))->default_value("vcs"),
       "name of the visible cross section graph (TGraphErrors*)")
      ("efficiency-name,e", po::value<std::string>(&(opts->efficiency_name)),
       "name of the detection efficiency object (TEfficiency*)")
      ( "ifname,i",
        po::value<std::string>(&(opts->ifname))->default_value("vcs.root"),
        "path to input file")
      ("ofname,o",
       po::value<std::string>(&(opts->ofname))->default_value("cond_num_test.root"),
       "path to output file")
      ("interp,r", po::value<std::string>(&(opts->interp)),
       "path to JSON file with interpolation settings");
}

void help(const po::options_description& desc) {
  std::cout << desc << std::endl;
}

double evalCondNumber(ISRSolverSLAE* solver) {
  solver->evalEqMatrix();
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(solver->getIntegralOperatorMatrix());
  double cond = svd.singularValues()(0) /
                svd.singularValues()(svd.singularValues().size()-1);
  return cond;
}

int main(int argc, char* argv[]) {
  po::options_description desc("This tool was designed to plot the condition number of an integral "
                               "operator matrix versus the center-of-mass energy spread. Allowed options");
  CmdOptions opts;
  setOptions(&desc, &opts);
  po::variables_map vmap;
  po::store(po::parse_command_line(argc, argv, desc), vmap);
  po::notify(vmap);
  if (vmap.count("help")) {
    help(desc);
    return 0;
  }
  ISRSolverSLAE solver(opts.ifname, {
      .efficiencyName = opts.efficiency_name,
      .visibleCSGraphName = opts.vcs_name,
      .thresholdEnergy = opts.thsd,
      .energyUnitMeVs = false});
  if (vmap.count("interp")) {
    solver.setRangeInterpSettings(opts.interp);
  }
  std::vector<double> ensigma;
  std::vector<double> condnums;
  ensigma.reserve(opts.n);
  condnums.reserve(opts.n);
  const double h = opts.ensigma_max / (opts.n - 1);
  ensigma.push_back(0);
  condnums.push_back(evalCondNumber(&solver));
  solver.enableEnergySpread();
  for (int i = 1; i < opts.n; ++i) {
    double sigma = i * h;
    ensigma.push_back(sigma);
    solver.resetECMErrors(sigma * Eigen::VectorXd::Ones(solver.ecm().size()));
    condnums.push_back(evalCondNumber(&solver));
  }
  TGraph g_cond(opts.n, ensigma.data(), condnums.data());
  auto fl = TFile::Open(opts.ofname.c_str(), "recreate");
  g_cond.Write("condnums");
  fl->Close();
  delete fl;
  return 0;
}
