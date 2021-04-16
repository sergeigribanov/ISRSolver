#include <iostream>
#include <algorithm>
#include <string>
#include <boost/program_options.hpp>
#include <TFile.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include "Interpolator.hpp"
namespace po = boost::program_options;

typedef struct {
  double thsd;
  std::string graph_name;
  std::string ifname;
  std::string ofname;
  std::string interp;
} CmdOptions;

void setOptions(po::options_description* desc, CmdOptions* opts) {
  desc->add_options()
      ("help,h", "Help.")
      ("thsd,t", po::value<double>(&(opts->thsd)), "Threshold (GeV).")
      ("cs-graph-name,c",
       po::value<std::string>(&(opts->graph_name))->default_value("bcs"),
       "Name of a cross section TGraphErrors.")
      ("ifname,i",
       po::value<std::string>(&(opts->ifname))->default_value("input.root"),
       "Path to input file.")
      ("ofname,o",
       po::value<std::string>(&(opts->ofname))->default_value("output.root"),
       "Path to output file.")
      ("interp,r",
       po::value<std::string>(&(opts->interp)),
       "Path to JSON file with interpolation settings.");
}

void help(const po::options_description& desc) {
  std::cout << desc << std::endl;
}

void interpBasis(const Eigen::VectorXd& y,
                 const Interpolator& interp,
                 double emax,
                 const CmdOptions& opts) {
  std::function<double(double*, double*)> fcn =
      [interp, y](double* px, double*) {
        double result = 0;
        for (int i = 0; i < y.rows(); ++i) {
          result += interp.basisEval(i, px[0]) * y(i);
        }
        return result;
      };
  TF1 f0("f_interp", fcn, opts.thsd, emax, 0);
  f0.SetNpx(10000);
  auto fl = TFile::Open(opts.ofname.c_str(), "recreate");
  fl->cd();
  f0.Write();
  fl->Close();
  delete fl;
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
  auto fl = TFile::Open(opts.ifname.c_str(), "read");
  auto bcs = dynamic_cast<TGraphErrors*>(fl->Get(opts.graph_name.c_str()));
  const int n = bcs->GetN();
  std::vector<double> x(n);
  std::vector<double> y(n);
  std::copy(bcs->GetX(), bcs->GetX() + n, x.begin());
  std::copy(bcs->GetY(), bcs->GetY() + n, y.begin());
  fl->Close();
  delete fl;
  Eigen::Map<Eigen::VectorXd> cmEnergies(x.data(), x.size());
  Eigen::Map<Eigen::VectorXd> v(y.data(), y.size());
  Interpolator* interp;
  if (vmap.count("interp")) {
    interp = new Interpolator(opts.interp, cmEnergies, opts.thsd);
  } else {
    interp = new Interpolator(cmEnergies, opts.thsd);
  }
  double emax = *std::max_element(x.begin(), x.end());
  interpBasis(v, *interp, emax, opts);
  delete interp;
  return 0;
}
