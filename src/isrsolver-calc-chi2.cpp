#include <iostream>
#include <utility>
#include <TFile.h>
#include <TF1.h>
#include <TMatrixD.h>
#include <TGraphErrors.h>
#include <Eigen/Dense>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

typedef struct {
  std::string name_of_model_bcs;
  std::string bcs_name;
  std::string path_to_model;
  std::string ifname;
} CmdOptions;

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

void help(const po::options_description& desc) {
  std::cout << desc << std::endl;
}

std::pair<double, double> evalChi2(const CmdOptions& opts) {
  auto ifl = TFile::Open(opts.ifname.c_str(), "read");
  ifl->cd();
  auto ibcs = dynamic_cast<TGraphErrors*>(ifl->Get(opts.bcs_name.c_str())->Clone("ibcs"));
  auto icovInv = dynamic_cast<TMatrixD*>(ifl->Get("invCovMatrixBornCS")->Clone("icovInv"));
  auto icov = dynamic_cast<TMatrixD*>(ifl->Get("covMatrixBornCS")->Clone("icov"));
  ifl->Close();

  int n = ibcs->GetN();
  Eigen::Map<Eigen::VectorXd> x(ibcs->GetX(), n);
  Eigen::Map<Eigen::VectorXd> data(ibcs->GetY(), n);
  Eigen::Map<Eigen::MatrixXd> invCov(icovInv->GetMatrixArray(), n, n);
  Eigen::Map<Eigen::MatrixXd> covMx(icov->GetMatrixArray(), n, n);
  Eigen::VectorXd mdata =Eigen::VectorXd::Zero(n);

  auto mfl = TFile::Open(opts.path_to_model.c_str(), "read");
  mfl->cd();
  auto mbcs = dynamic_cast<TF1*>(mfl->Get(opts.name_of_model_bcs.c_str()));
  for (int i = 0; i < x.size(); ++i) {
    mdata(i) = mbcs->Eval(x(i));
  }
  mfl->Close();
  delete mfl;
  Eigen::VectorXd dy = data - mdata;
  double chi2 = dy.dot(invCov * dy);
  Eigen::VectorXd diagErr2 = covMx.diagonal().array().pow(-1.);
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
