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

/**
 * Chi-square test
 * @param solver a solver
 * @param args a chi-square test arguments
 * @param vcs a visible cross section
 * @param bcs0 a reference Born cross section
 * @param vcsErr a visible cross section errors
 * @see ISRSolverSLE
 * @see Chi2TestArgs
 */
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
    /**
     * Random redraw of the initial visible cross section
     */
    solver->resetVisibleCS(randomDrawVisCS(vcs, vcsErr));
    /**
     * Find a numerical solution
     */
    solver->solve();
    /**
     * Find the difference between the numerical solution and
     the model (or initial) Born cross section
     */
    dbcs = solver->bcs() - bcs0;
    /**
     * Calculating chi-square
     */
    double tmp_chi2 = dbcs.dot(lu.solve(dbcs));
    /**
     * Push the chi-square value to a temporary vector
     */
    chi2s.push_back(tmp_chi2);
  }
  /**
   * Getting the lowest chi-square value
   */
  const double chi2Min = *std::min_element(chi2s.begin(), chi2s.end());
  /**
   * Getting the highest chi-square value
   */
  const double chi2Max = *std::max_element(chi2s.begin(), chi2s.end());
  /**
   * Creating chi-square histogram
   */
  TH1F chi2Hist("chi2Hist", "", 128, chi2Min, chi2Max);
  /**
   * Fill chi-square histogram
   */
  for (const double value : chi2s) {
    chi2Hist.Fill(value);
  }
  /**
   * Creating chi-square fit function
   */
  std::function<double(double*, double*)> fcn_chi2 =
      [](double* px, double* par) {
        return par[0] * ROOT::Math::chisquared_pdf(px[0], par[1]);
      };
  /**
   * Creating chi-square fit function in a form of TF1 object
   */
  TF1 f1("f1", fcn_chi2, chi2Min, chi2Max, 2);
  f1.SetParameter(0, args.initialChi2Ampl);
  f1.SetParameter(1, vcs.size());
  f1.SetNpx(10000);
  /**
   * Fitting to the chi-square histogram
   */
  chi2Hist.Fit(&f1, "LME+");
  /**
   * Saving chi-square histogram and fit function in a file
   */
  auto fl = TFile::Open(args.outputPath.c_str(), "recreate");
  chi2Hist.Write();
  fl->Close();
  delete fl;
}

/**
 * Chi-square model test
 * @parm solver a solver
 * @parm modelArgs a structure with input arguments
 * @see ISRSolverSLE
 * @see Chi2TestModelArgs
 */
void chi2TestModel(ISRSolverSLE* solver,
                   const Chi2TestModelArgs& args) {
  /**
   * Finding a numerical solution
   */
  solver->solve();
  Eigen::VectorXd vcsErr = solver->vcsErr();
  Eigen::VectorXd vcs = Eigen::VectorXd(vcsErr.size());
  Eigen::VectorXd bcs0 = Eigen::VectorXd(vcsErr.size());
  /**
   * Reading a file with model data
   */
  auto ifl = TFile::Open(args.modelPath.c_str(), "read");
  /**
   * Getting a pointer to a model visible cross section graph
   */
  auto g_bcs = dynamic_cast<TGraphErrors*>(ifl->Get(args.modelBCSName.c_str()));
  /**
   * Getting a pointer to a model Born cross section graph
   */
  auto g_vcs = dynamic_cast<TGraphErrors*>(ifl->Get(args.modelVCSName.c_str()));
  // !!! TO DO : exception if vcs and bcs0 sizes are different
  // TO DO : sort array, copy buffers
  /**
   * Fill model visible and Born cross section vectors
   */
  // Eigen::Map<Eigen::VectorXd> vcs(g_vcs->GetY(), g_vcs->GetN());
  // Eigen::Map<Eigen::VectorXd> bcs0(g_bcs->GetY(), g_bcs->GetN());
  for (int i = 0; i < g_bcs->GetN(); ++i) {
    bcs0(i) = g_bcs->GetY()[i];
    vcs(i) = g_vcs->GetY()[i];
  }
  ifl->Close();
  delete ifl;
  /**
   * Run chi-square test
   */
  chi2Test(solver,
           {.n = args.n,
            .initialChi2Ampl = args.initialChi2Ampl,
            .outputPath = args.outputPath},
           vcs, bcs0, vcsErr);
}

/**
 * Chi-square data test
 * @param solver a solver
 * @param args a structure with input arguments
 * @see ISRSolverSLE
 * @see Chi2TestArgs
 */
void chi2TestData(ISRSolverSLE* solver,
                  const Chi2TestArgs& args) {
  /**
   * Finding a numerical solution
   */
  solver->solve();
  Eigen::VectorXd vcsErr = solver->vcsErr();
  Eigen::VectorXd vcs = Eigen::VectorXd(vcsErr.size());
  Eigen::VectorXd bcs0 = Eigen::VectorXd(vcsErr.size());
  vcs = solver->vcs();
  /**
   * Accumulator used to obtain a mean numerical solution
   */
  std::vector<bacc::accumulator_set<double, bacc::stats<bacc::tag::mean>>> accs(bcs0.size());
  /**
   * Obtaining a mean numerical solution
   */
  for (int i = 0; i < args.n; ++i) {
    solver->resetVisibleCS(randomDrawVisCS(vcs, vcsErr));
    solver->solve();
    for (int j = 0; j < bcs0.size(); ++j) {
      accs[j](solver->bcs()(j));
    }
  }
  for (int i = 0; i < bcs0.size(); ++i) {
    /**
     * Populating the mean values of a numerical solution
     */
    bcs0(i) = bacc::mean(accs[i]);
  }
  /**
   * Run chi-square test
   */
  chi2Test(solver, args,
           vcs, bcs0, vcsErr);
}
