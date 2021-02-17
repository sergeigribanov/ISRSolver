#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <functional>
#include <Eigen/Core>
#include <Eigen/SVD>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TFile.h>
#include <TRandom3.h>
#include <Math/PdfFuncMathCore.h>
#include "ISRSolverSLAE.hpp"
#include "ISRSolverTikhonov.hpp"
#include "ISRSolverTSVD.hpp"
#include "utils.hpp"
namespace po = boost::program_options;

typedef struct {
  double thsd;
  double alpha;
  int k;
  int n;
  std::string path_to_model;
  std::string name_of_model_bcs;
  std::string name_of_model_vcs;
  std::string vcs_name;
  std::string efficiency_name;
  std::string ifname;
  std::string ofname;
  std::string interp;
  std::string solver;
} CmdOptions;

void setOptions(po::options_description* desc, CmdOptions* opts) {
  desc->add_options()("help,h",
                      "A simple tool designed in order to find numerical"
                      "solution of the Kuraev-Fadin equation.")
      ("thsd,t", po::value<double>(&(opts->thsd)), "Threshold (GeV).")
      ("num-rnd-draws,n", po::value<int>(&(opts->n))->default_value(100),
       "Number of visible cross section random draws.")
      ("upper-tsvd-index,k", po::value<int>(&(opts->k)), "Upper TSVD index")
      ("keep-one,z", "Keep only k-th SVD harmonic")
      ("enable-energy-spread,g", "Enable energy spread")
      // ("enable-solution-positivity,p", "Setting positive limits to solution")
      ("disable-solution-norm,f", 
       "Disable solution norm in Tihonov's regularization functional")
      ("disable-solution-derivative-norm,d",
       "Disable solution derivative norm in Tikhonov's regularization functional")
      ("alpha,a", po::value<double>(&(opts->alpha)),
       "Thikhonov's regularization parameter.")(
           "solver,s", po::value<std::string>(&(opts->solver)),
           "Solver: SLAE, Tikhonov")(
               "vcs-name,v", po::value<std::string>(&(opts->vcs_name))->default_value("vcs"),
               "Name of the visible cross section graph.")
      ("efficiency-name,e", po::value<std::string>(&(opts->efficiency_name)),
       "TEfficiency object name")
      ("use-model,u", po::value<std::string>(&(opts->path_to_model)),
       "Path to the file with the model Born and visible cross section TGraphErrors (if needed)")
      ("model-bcs-name,b", po::value<std::string>(&(opts->name_of_model_bcs))->default_value("bcsSRC"),
       "Name of the model Born cross section TGraphErrors function")
      ("model-vcs-name,c", po::value<std::string>(&(opts->name_of_model_vcs))->default_value("vcsBlured"),
       "Name of the model visible cross section TGraphErrors function")
      ( "ifname,i",
        po::value<std::string>(&(opts->ifname))->default_value("vcs.root"),
        "Path to input file.")(
            "ofname,o",
            po::value<std::string>(&(opts->ofname))->default_value("isrsolver-chi2-distribution.root"),
            "Path to output file.")("interp,r",
                                    po::value<std::string>(&(opts->interp)),
                                    "Path to JSON file with interpolation settings.");
}

void help(const po::options_description& desc) {
  std::cout << desc << std::endl;
}

Eigen::VectorXd randomDrawVisCS(const Eigen::VectorXd& vcs,
                                const Eigen::VectorXd& vcsErr) {
  Eigen::VectorXd result = vcs;
  for (int i = 0; i < result.size(); ++i) {
    result(i) = gRandom->Gaus(result(i), vcsErr(i));
  }
  return result;
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
  if (!vmap.count("use-model")) {
    std::cout << "[!] use-model command option is obligatory" << std::endl;
    return 0;
  }
  BaseISRSolver* solver = nullptr;
  if (opts.solver == "SLAE") {
    solver = new ISRSolverSLAE(opts.ifname, {
	.efficiencyName = opts.efficiency_name,
	.visibleCSGraphName = opts.vcs_name,
	.thresholdEnergy = opts.thsd,
	.energyUnitMeVs = false});
  } else if (opts.solver == "Tikhonov") {
    solver =
        new ISRSolverTikhonov(opts.ifname, {
	    .efficiencyName = opts.efficiency_name,
	    .visibleCSGraphName = opts.vcs_name,
	    .thresholdEnergy = opts.thsd,
	    .energyUnitMeVs = false});
  } else if (opts.solver == "TSVD") {
    solver =
        new ISRSolverTSVD(opts.ifname, {
            .efficiencyName = opts.efficiency_name,
            .visibleCSGraphName = opts.vcs_name,
            .thresholdEnergy = opts.thsd,
            .energyUnitMeVs = false}, 1);
  }
  if (!solver) {
    std::cerr << "[!] Solver is not set." << std::endl;
    return 1;
  }
  ISRSolverSLAE* solverSLAE = dynamic_cast<ISRSolverSLAE*>(solver);
  if (!solverSLAE) {
    std::cout << "Wrong solver type!" << std::endl;
    return 0;
  }
  if (vmap.count("interp") && solverSLAE) {
    solverSLAE->setRangeInterpSettings(opts.interp);
  }
  if (vmap.count("enable-energy-spread")) {
    solverSLAE->enableEnergySpread();
  }
  ISRSolverTikhonov* solverTikhonov = dynamic_cast<ISRSolverTikhonov*>(solver);
  if (vmap.count("alpha") && solverTikhonov) {
    solverTikhonov->setAlpha(opts.alpha);
  }
  ISRSolverTSVD* solverTSVD = dynamic_cast<ISRSolverTSVD*>(solver);
  if (vmap.count("upper-tsvd-index") && solverTSVD) {
    solverTSVD->setTruncIndexUpperLimit(opts.k);
  }
  if (vmap.count("keep-one") && solverTSVD) {
    solverTSVD->enableKeepOne();
  }
  if (vmap.count("disable-solution-norm") && solverTikhonov) {
    solverTikhonov->disableSolutionNorm2();
  }
  if (vmap.count("disable-solution-derivative-norm") && solverTikhonov) {
    solverTikhonov->disableSolutionDerivativeNorm2();
  }
  solver->solve();
  Eigen::VectorXd ecm = solver->ecm();
  Eigen::VectorXd vcsErr = solver->vcsErr();
  Eigen::VectorXd vcs = Eigen::VectorXd(vcsErr.size());
  Eigen::VectorXd bcs0 = Eigen::VectorXd(vcsErr.size());
  auto ifl = TFile::Open(opts.path_to_model.c_str(), "read");
  auto g_bcs = dynamic_cast<TGraphErrors*>(ifl->Get(opts.name_of_model_bcs.c_str()));
  auto g_vcs = dynamic_cast<TGraphErrors*>(ifl->Get(opts.name_of_model_vcs.c_str()));
  // TO DO : sort array, copy buffers
  for (int i = 0; i < g_bcs->GetN(); ++i) {
    bcs0(i) = g_bcs->GetY()[i];
    vcs(i) = g_vcs->GetY()[i];
  }
  ifl->Close();
  delete ifl;
  Eigen::VectorXd dbcs;
  std::vector<std::shared_ptr<TH1F>> bias;
  bias.reserve(bcs0.size());
  for (int i = 0; i < bcs0.size(); ++i) {
    std::string hist_name = boost::str(boost::format("bias_pt_%1%") % i);
    bias.push_back(std::shared_ptr<TH1F>(
        new TH1F(hist_name.c_str(), "", 1024, -2, 2)));
  }
  for (int i = 0; i < opts.n; ++i) {
    solver->resetVisibleCS(randomDrawVisCS(vcs, vcsErr));
    solver->solve();
    dbcs = solver->bcs() - bcs0;
    for (int j = 0; j < dbcs.size(); ++j) {
      bias[j].get()->Fill(dbcs(j) / bcs0(j));
    }
  }
  std::vector<double> mean_values;
  std::vector<double> mean_errors;
  mean_values.reserve(ecm.size());
  mean_errors.reserve(ecm.size());
  for (int i = 0; i < ecm.size(); ++i) {
    mean_values.push_back(bias[i].get()->GetMean());
    mean_errors.push_back(bias[i].get()->GetMeanError());
  }
  TGraphErrors means(ecm.size(), ecm.data(), mean_values.data(), 0, mean_errors.data());
  auto fl = TFile::Open(opts.ofname.c_str(), "recreate");
  means.Write("means");
  for (int i = 0; i < bcs0.size(); ++i) {
    bias[i].get()->Write();
  }
  fl->Close();
  delete solver;
  return 0;
}
