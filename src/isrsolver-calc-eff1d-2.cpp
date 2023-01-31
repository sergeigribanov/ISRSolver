#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <functional>
#include <boost/program_options.hpp>
#include <nlohmann/json.hpp>
#include <TFile.h>
#include <TF1.h>
#include <TGraph.h>
#include <TH1D.h>
#include "KuraevFadin.hpp"
namespace po = boost::program_options;
using json = nlohmann::json;

/**
 * A part of program options
 */
typedef struct {
  /**
   * Threshold energy
   */
  double thsd;
  /**
   * Name of the Born cross section graph
   * (TF1)
   */
  std::string bcs_name;
  /**
   * Name of the detection efficiency object
   * (TH1D)
   */
  std::string efficiency_name;
  std::string eff_json_path;
  std::string en_json_path;
  /**
   * Path to that input .root file that contains
   * the Born cross section
   */
  std::string ifname;
  /**
   * Output file path
   */
  std::string ofname;
} CmdOptions;

/**
 * Setting up program options
 */
void setOptions(po::options_description* desc, CmdOptions* opts) {
  desc->add_options()
      ("help,h", "help message")
      ("thsd,t", po::value<double>(&(opts->thsd)), "threshold energy (GeV)")
      ("bcs-name,b", po::value<std::string>(&(opts->bcs_name))->default_value("f_bcs"),
       "name of a Born cross section graph (TF1)")
      ( "ifname,i",
        po::value<std::string>(&(opts->ifname))->default_value("bcs.root"),
        "path to input file")
      ("ofname,o",
       po::value<std::string>(&(opts->ofname))->default_value("bcs.root"),
       "path to output file")
      ("efficiency-name,e", po::value<std::string>(&(opts->efficiency_name))->default_value("eps_x"),
       "name of a detection efficiency histogram")
      ("eff-json-path,j", po::value<std::string>(&(opts->eff_json_path)),
       "path to json-file with efficiency file-names")
      ("en-json-path,p", po::value<std::string>(&(opts->en_json_path)),
      "path to json-file with energy points");
}

/**
 * Help message
 */
void help(const po::options_description& desc) {
  std::cout << desc << std::endl;
}

int main(int argc, char* argv[]) {
  po::options_description desc("   Solver that uses the naive method. The naive method consists in reducing the"
                               "integral equation to a system of linear differential equations. Allowed options");
  CmdOptions opts;
  setOptions(&desc, &opts);
  po::variables_map vmap;
  po::store(po::parse_command_line(argc, argv, desc), vmap);
  po::notify(vmap);
  if (vmap.count("help")) {
    help(desc);
    return 0;
  }
  std::vector<std::string> hist_paths;
  if (vmap.count("eff-json-path")) {
    std::ifstream f_json_eff(opts.eff_json_path);
    json eff_data = json::parse(f_json_eff);
    for (std::size_t i = 0; i < eff_data.size(); ++i) {
      hist_paths.push_back(eff_data[i]);
    }
  }
  std::vector<double> ens;
  if (vmap.count("en-json-path")) {
    std::ifstream f_json_en(opts.en_json_path);
    json en_data = json::parse(f_json_en);
    for (std::size_t i = 0; i < en_data.size(); ++i) {
      ens.push_back(en_data[i]["energy"]);
    }
  }
  assert(ens.size() == hist_paths.size());
  std::vector<double> effs;
  auto ifl = TFile::Open(opts.ifname.c_str(), "read");
  auto tf1_bcs = dynamic_cast<TF1*>(ifl->Get(opts.bcs_name.c_str())->Clone());
  ifl->Close();
  delete ifl;
  std::function<double(double)> f_bcs = [opts, tf1_bcs](double en) {
    if (en <= opts.thsd) {
      return 0.;
    }
    return tf1_bcs->Eval(en);
  };
  const double s_t = opts.thsd * opts.thsd;
  for (std::size_t i = 0; i < ens.size(); ++i) {
    auto efl = TFile::Open(hist_paths[i].c_str(), "read");
    auto eff_2d = dynamic_cast<TH1D*>(efl->Get(opts.efficiency_name.c_str())->Clone());
    eff_2d->SetDirectory(0);
    efl->Close();
    delete efl;
    std::function<double(double)> eff_1d = [eff_2d](double x) {
      int bin = eff_2d->FindBin(x);
      return eff_2d->GetBinContent(bin);
    };
    const double s = ens[i] * ens[i];
    const double conv0 = convolutionKuraevFadin(ens[i], f_bcs, 0., 1. - s_t / s);
    const double conv1 = convolutionKuraevFadin_1d(ens[i], f_bcs, 0., 1. - s_t / s, eff_1d);
    effs.push_back(conv1 / conv0);
    delete eff_2d;
  }
  auto ofl = TFile::Open(opts.ofname.c_str(), "recreate");
  TGraph gr(ens.size(), ens.data(), effs.data());
  gr.Write("eff_1d_calc");
  ofl->Close();
  delete ofl;
  return 0;
}
