#include <vector>
#include <algorithm>
#include <functional>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <Eigen/Core>
#include <Eigen/SVD>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TRandom3.h>
#include <Math/PdfFuncMathCore.h>
#include "Chi2Test.hpp"
#include "Utils.hpp"
namespace bacc = boost::accumulators;

void chi2Test(ISRSolverSLE* solver,
              const Chi2TestArgs& args,
              const Eigen::VectorXd& vcs,
              const Eigen::VectorXd& bcs0,
              const Eigen::VectorXd& vcsErr) {
  Eigen::VectorXd dbcs;
  std::vector<double> chi2s;
  chi2s.reserve(args.n);
  Eigen::FullPivLU<Eigen::MatrixXd> lu(solver->getBornCSCovMatrix());
  for (int i = 0; i < args.n; ++i) {
    solver->resetVisibleCS(randomDrawVisCS(vcs, vcsErr));
    solver->solve();
    dbcs = solver->bcs() - bcs0;
    double tmp_chi2 = dbcs.dot(lu.solve(dbcs));
    chi2s.push_back(tmp_chi2);
  }
  const double chi2Min = *std::min_element(chi2s.begin(), chi2s.end());
  const double chi2Max = *std::max_element(chi2s.begin(), chi2s.end());
  TH1F chi2Hist("chi2Hist", "", 128, chi2Min, chi2Max);
  for (const double value : chi2s) {
    chi2Hist.Fill(value);
  }
  std::function<double(double*, double*)> fcn_chi2 =
      [](double* px, double* par) {
        return par[0] * ROOT::Math::chisquared_pdf(px[0], par[1]);
      };
  TF1 f1("f1", fcn_chi2, chi2Min, chi2Max, 2);
  f1.SetParameter(0, args.initialChi2Ampl);
  f1.SetParameter(1, vcs.size());
  f1.SetNpx(10000);
  chi2Hist.Fit(&f1, "LME+");
  auto fl = TFile::Open(args.outputPath.c_str(), "recreate");
  chi2Hist.Write();
  fl->Close();
  delete fl;
}

void chi2TestModel(ISRSolverSLE* solver,
                   const Chi2TestModelArgs& args) {
  solver->solve();
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
  chi2Test(solver,
           {.n = args.n,
            .initialChi2Ampl = args.initialChi2Ampl,
            .outputPath = args.outputPath},
           vcs, bcs0, vcsErr);
}


void chi2TestData(ISRSolverSLE* solver,
                  const Chi2TestArgs& args) {
  solver->solve();
  Eigen::VectorXd vcsErr = solver->vcsErr();
  Eigen::VectorXd vcs = Eigen::VectorXd(vcsErr.size());
  Eigen::VectorXd bcs0 = Eigen::VectorXd(vcsErr.size());
  vcs = solver->vcs();
  std::vector<bacc::accumulator_set<double, bacc::stats<bacc::tag::mean>>> accs(bcs0.size());
  for (int i = 0; i < args.n; ++i) {
    solver->resetVisibleCS(randomDrawVisCS(vcs, vcsErr));
    solver->solve();
    for (int j = 0; j < bcs0.size(); ++j) {
      accs[j](solver->bcs()(j));
    }
  }
  for (int i = 0; i < bcs0.size(); ++i) {
    bcs0(i) = bacc::mean(accs[i]);
  }
  chi2Test(solver, args,
           vcs, bcs0, vcsErr);
}
