#include <memory>
#include <vector>
#include <boost/format.hpp>
#include <Eigen/Dense>
#include <TFile.h>
#include <TF1.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TRandom3.h>
#include "RatioTest.hpp"
#include "Utils.hpp"

void ratioTestModel(ISRSolverSLE* solver,
                    const RatioTestModelArgs& args) {
  solver->solve();
  Eigen::VectorXd ecm = solver->ecm();
  Eigen::VectorXd vcsErr = solver->vcsErr();
  Eigen::VectorXd vcs = Eigen::VectorXd(vcsErr.size());
  Eigen::VectorXd bcs0 = Eigen::VectorXd(vcsErr.size());
  auto ifl = TFile::Open(args.modelPath.c_str(), "read");
  auto g_bcs = dynamic_cast<TGraphErrors*>(ifl->Get(args.modelBCSName.c_str()));
  auto g_vcs = dynamic_cast<TGraphErrors*>(ifl->Get(args.modelVCSName.c_str()));
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
  for (int i = 0; i < args.n; ++i) {
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
  auto fl = TFile::Open(args.outputPath.c_str(), "recreate");
  means.Write("means");
  for (int i = 0; i < bcs0.size(); ++i) {
    bias[i].get()->Write();
  }
  fl->Close();
  delete fl;
}
