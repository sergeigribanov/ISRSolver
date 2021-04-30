#include <iostream>
#include <utility>
#include <boost/program_options.hpp>
#include <Eigen/Dense>
#include <TFile.h>
#include <TF1.h>
#include <TMatrixD.h>
#include <TGraphErrors.h>
namespace po = boost::program_options;

/**
 * A part of program options
 */
typedef struct {
  /**
   * Name of model Born cross section object (TF1)
   */
  std::string name_of_model_bcs;
  /**
   * Name of the numerical solution (Born cross section)
   * object (TGraphErrors)
   */
  std::string bcs_name;
  /**
   * Path to a .root file with a model data
   */
  std::string path_to_model;
  /**
   * Path to a .root file with a numerical
   * solution (Born cross section) data
   */
  std::string ifname;
} CmdOptions;

/**
 * Setting up program options
 */
void setOptions(po::options_description* desc, CmdOptions* opts) {
  desc->add_options()
      ("help,h", "help message")
      ("bcs-name,b",
       po::value<std::string>(&(opts->bcs_name))->default_value("bcs"),
       "name of numerical solution (Born cross section)")
      ("use-model,u", po::value<std::string>(&(opts->path_to_model)),
       "path to the file with the model Born cross section graph (TGraphErrors*)")
      ("model-bcs-name,m",
       po::value<std::string>(&(opts->name_of_model_bcs))->default_value("f_bcs"),
       "name of the model Born cross section function (TF1*)")
      ( "ifname,i",
        po::value<std::string>(&(opts->ifname))->default_value("bcs.root"),
        "path to input file");
}

/**
 * Help message
 */
void help(const po::options_description& desc) {
  std::cout << desc << std::endl;
}

/**
 * Evaluate chi-square
 * @param opts a part of the program options
 */
std::pair<double, double> evalChi2(const CmdOptions& opts) {
  /**
   * Opening input file that contains Born cross section (numerical solution),
   * its covariance matrix and inverse covariance matrix
   */
  auto ifl = TFile::Open(opts.ifname.c_str(), "read");
  ifl->cd();
  /**
   * Reading Born cross section (numerical solution)
   */
  auto ibcs = dynamic_cast<TGraphErrors*>(ifl->Get(opts.bcs_name.c_str())->Clone("ibcs"));
  /**
   * Reading inverse covariance matrix of the numerical solution
   */
  auto icovInv = dynamic_cast<TMatrixD*>(ifl->Get("invCovMatrixBornCS")->Clone("icovInv"));
  /**
   * Reading covariance matrix of the numerical solution
   */
  auto icov = dynamic_cast<TMatrixD*>(ifl->Get("covMatrixBornCS")->Clone("icov"));
  ifl->Close();

  int n = ibcs->GetN();
  Eigen::Map<Eigen::VectorXd> x(ibcs->GetX(), n);
  Eigen::Map<Eigen::VectorXd> data(ibcs->GetY(), n);
  Eigen::Map<Eigen::MatrixXd> invCov(icovInv->GetMatrixArray(), n, n);
  Eigen::Map<Eigen::MatrixXd> covMx(icov->GetMatrixArray(), n, n);
  Eigen::VectorXd mdata =Eigen::VectorXd::Zero(n);

  /**
   * Opening the file with a model data
   */
  auto mfl = TFile::Open(opts.path_to_model.c_str(), "read");
  mfl->cd();
  /**
   * Reading model Born cross section
   */
  auto mbcs = dynamic_cast<TF1*>(mfl->Get(opts.name_of_model_bcs.c_str()));
  /**
   * Converting model Born cross section to a vector format
   */
  for (int i = 0; i < x.size(); ++i) {
    mdata(i) = mbcs->Eval(x(i));
  }
  mfl->Close();
  delete mfl;
  Eigen::VectorXd dy = data - mdata;
  /**
   * Calculating chi-square using full covariance matrix
   */
  double chi2 = dy.dot(invCov * dy);
  Eigen::VectorXd diagErr2 = covMx.diagonal().array().pow(-1.);
  /**
   * Calculating chi-square using only diagonal elements of
   * covariance matrix
   */
  double chi2Diag = dy.dot(diagErr2.asDiagonal() * dy);
  return std::make_pair(chi2, chi2Diag);
}

int main(int argc, char* argv[]) {
  po::options_description desc("   Chi-square calculation");
  CmdOptions opts;
  setOptions(&desc, &opts);
  po::variables_map vmap;
  po::store(po::parse_command_line(argc, argv, desc), vmap);
  po::notify(vmap);
  if (vmap.count("help")) {
    help(desc);
    return 0;
  }
  auto chi2s = evalChi2(opts);
  std::cout << "chi-square = " << chi2s.first << std::endl;
  std::cout << "diagonal chi-square = " << chi2s.second << std::endl;
  return 0;
}
