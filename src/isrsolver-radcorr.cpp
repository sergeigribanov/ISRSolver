#include <iostream>
#include <vector>
#include <functional>
#include <boost/program_options.hpp>
#include <TFile.h>
#include <TF1.h>
#include <TGraph.h>
#include "kuraev_fadin.hpp"
namespace po = boost::program_options;

typedef struct {
  std::size_t n;
  double thsd;
  double minen;
  double maxen;
  std::string ifname;
  std::string bcs_fcn_name;
  std::string ofname;
} CmdOptions;

void setOptions(po::options_description* desc, CmdOptions* opts) {
  desc->add_options()
      ("help,h", "Radiative correction calculator.")
      ("thsd,t", po::value<double>(&(opts->thsd)), "Threshold (GeV).")
      ("number-of-radocrr-points,n", po::value<std::size_t>(&(opts->n))->default_value(1000),
       "Number of radiative correction points.")
      ("minen,m",  po::value<double>(&(opts->minen)), "Minimum energy.")
      ("maxen,x",  po::value<double>(&(opts->maxen)), "Maximum energy.")
      ( "ifname,i",
        po::value<std::string>(&(opts->ifname))->default_value("vcs.root"),
        "Path to input file.")
      ("ofname,o",
       po::value<std::string>(&(opts->ofname))->default_value("bcs.root"),
       "Path to output file.")
      ("born-cs-fcn-name",
       po::value<std::string>(&(opts->bcs_fcn_name))->default_value("f_bcs"),
       "The name of the Born cross section TF1 function.");
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
  if (!vmap.count("maxen") ||
      !vmap.count("minen") ||
      !vmap.count("thsd")) {
    std::cout << "[!] You need to set threshold energy, maximum and minimum energy" << std::endl;
    return 0;
  }
  auto fl = TFile::Open(opts.ifname.c_str(), "read");
  auto fbcs = dynamic_cast<TF1*>(fl->Get(opts.bcs_fcn_name.c_str())->Clone());
  fl->Close();
  delete fl;
  std::function<double(double)> fcn =
      [fbcs](double energy) {
        double result = fbcs->Eval(energy);
        return result;
      };

  std::vector<double> ens;
  std::vector<double> radcorrs;
  ens.reserve(opts.n);
  radcorrs.reserve(opts.n);
  const double s_th = opts.thsd * opts.thsd;
  const double eh = (opts.maxen - opts.minen) / opts.n;
  for (std::size_t i = 1; i < opts.n + 1; ++i) {
    double en = opts.minen + i * eh;
    ens.push_back(en);
    radcorrs.push_back(kuraev_fadin_convolution(en, fcn, 0, 1 - s_th / en / en) / fcn(en) - 1);
  }
  delete fbcs;
  TGraph gradcorr(opts.n, ens.data(), radcorrs.data());
  auto ofl = TFile::Open(opts.ofname.c_str(), "recreate");
  ofl->cd();
  gradcorr.Write("radcorr");
  ofl->Close();
  delete ofl;
  return 0;
}