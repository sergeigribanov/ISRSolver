#include <TF1.h>
#include <TFile.h>

#include <algorithm>
#include <boost/program_options.hpp>
#include <functional>
#include <iostream>
#include <string>

#include "kuraev_fadin.hpp"
#include "utils.hpp"
namespace po = boost::program_options;

typedef struct {
  double thsd;
  std::string fcn;
  std::string ifname;
  std::string ofname;
} CmdOptions;

void setOptions(po::options_description* desc, CmdOptions* opts) {
  desc->add_options()("help",
                      "A simple tool designed in order to find convolution of "
                      "a Born cross section "
                      "with the Kuraev-Fadin kernel.")(
      "thsd", po::value<double>(&(opts->thsd)), "Threshold (GeV).")(
      "fcn", po::value<std::string>(&(opts->fcn)),
      "fcn is a name of a function used to describe a Born cross section, "
      "which will be convaluted with the Kuraev-Fadin kernel.")(
      "ifname",
      po::value<std::string>(&(opts->ifname))->default_value("input.root"),
      "Path to input file.")(
      "ofname",
      po::value<std::string>(&(opts->ofname))->default_value("output.root"),
      "Path to output file.");
}

void help(const po::options_description& desc) {
  std::cout << desc << std::endl;
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
  auto fcn = find_object<TF1>(fl, opts.fcn);
  std::function<double(double)> born_fcn = [&opts, &fcn](double s) {
    double e0 = opts.thsd;
    double e1 = fcn->GetXmin();
    double en = std::sqrt(s);
    if (e0 < e1 && en < e1) {
      return fcn->Eval(e1) * (en - e0) / (e1 - e0);
    }
    return fcn->Eval(en);
  };
  std::function<double(double*, double*)> vcs_fcn =
      [&opts, &born_fcn](double* x, double* par) {
        double s = x[0] * x[0];
        double s_threshold = opts.thsd * opts.thsd;
        return kuraev_fadin_convolution(s, born_fcn, 0, 1 - s_threshold / s);
      };
  auto f_vcs = new TF1("f_vcs", vcs_fcn, std::max(opts.thsd, fcn->GetXmin()),
                       fcn->GetXmax(), 0);
  fl->Close();
  delete fl;
  auto fl_out = TFile::Open(opts.ofname.c_str(), "recreate");
  fl_out->cd();
  f_vcs->Write();
  fl_out->Close();
  delete f_vcs;
  delete fl_out;
  return 0;
}
