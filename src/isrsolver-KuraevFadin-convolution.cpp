#include <algorithm>
#include <functional>
#include <iostream>
#include <string>
#include <boost/program_options.hpp>
#include <TF1.h>
#include <TFile.h>
#include <TEfficiency.h>
#include "KuraevFadin.hpp"
#include "Utils.hpp"
namespace po = boost::program_options;

//!!! TO DO: insert efficiency
/**
 * A part of program options
 */
typedef struct {
  /**
   * Number of points used to plot the function
   */
  int n;
  /**
   * Threshold energy
   */
  double thsd;
  /**
   * Name of a function (TF1) to be convoluted with kernel function
   */
  std::string fcn;
  /**
   * Name of a detection efficiency object (TEfficiency)
   */
  std::string efficiency_name;
  /**
   * Path to the input .root file that contains function to be convoluted and
   * detection efficiency
   */
  std::string ifname;
  /**
   * Path to the output .root file
   */
  std::string ofname;
} CmdOptions;

/**
 * Setting up program options
 */
void setOptions(po::options_description* desc, CmdOptions* opts) {
  desc->add_options()("help,h", "help message")
      ("num-of-points,n", po::value<int>(&(opts->n))->default_value(10000), "number of center-of-mass energy points")
      ("thsd,t", po::value<double>(&(opts->thsd)), "threshold energy (GeV)")(
      "fcn,f", po::value<std::string>(&(opts->fcn)),
      "name of the function that to be convoluted with the kernel function F(x,s)")(
      "ifname,i",
      po::value<std::string>(&(opts->ifname))->default_value("input.root"),
      "path to input file")
      ("ofname,o",
      po::value<std::string>(&(opts->ofname))->default_value("output.root"),
      "path to output file")
      ("efficiency-name,e", po::value<std::string>(&(opts->efficiency_name)),
       "name of a detection efficiency object (TEfficiency*)");
}

/**
 * Help message
 */
void help(const po::options_description& desc) {
  std::cout << desc << std::endl;
}

int main(int argc, char* argv[]) {
  po::options_description desc(
      "  This tool is designed in order to make convolution of a some "
      "custom (TF1*) function with the kernel function F(x,s). Allowed options");
  CmdOptions opts;
  setOptions(&desc, &opts);
  po::variables_map vmap;
  po::store(po::parse_command_line(argc, argv, desc), vmap);
  po::notify(vmap);
  if (vmap.count("help")) {
    help(desc);
    return 0;
  }
  TEfficiency* teff = nullptr;
  /**
   * Opening the input file
   */
  auto fl = TFile::Open(opts.ifname.c_str(), "read");
  /**
   * Reading the function
   */
  auto fcn = dynamic_cast<TF1*>(fl->Get(opts.fcn.c_str()));
  /**
   * Reading the detection efficiency
   */
  if (vmap.count("efficiency-name")) {
    teff = dynamic_cast<TEfficiency*>(fl->Get(opts.efficiency_name.c_str())->Clone());
  }
  /**
   * Converting the function to a form of std::function
   */
  std::function<double(double)> born_fcn = [&opts, &fcn](double s) {
    double e0 = opts.thsd;
    double e1 = fcn->GetXmin();
    double en = std::sqrt(s);
    if (e0 < e1 && en < e1) {
      return fcn->Eval(e1) * (en - e0) / (e1 - e0);
    }
    return fcn->Eval(en);
  };
  /**
   * Converting the detection efficiency to a form of std::function
   */
  std::function<double(double, double)> eff =
      [teff](double x, double en) {
        int bin = teff->FindFixBin(x, en);
        double result = teff->GetEfficiency(bin);
        return result;
      };
  /**
   * Creating convolution function
   */
  std::function<double(double*, double*)> vcs_fcn =
      [opts, born_fcn, teff, eff](double* x, double* par) {
        double s = x[0] * x[0];
        double s_threshold = opts.thsd * opts.thsd;
        double result = 0;
        if (teff) {
          result = convolutionKuraevFadin(s, born_fcn, 0, 1 - s_threshold / s, eff);
        } else {
          result = convolutionKuraevFadin(s, born_fcn, 0, 1 - s_threshold / s);
        }
        return result;
      };
  auto f_vcs = new TF1("f_vcs", vcs_fcn, std::max(opts.thsd, fcn->GetXmin()),
                       fcn->GetXmax(), 0);
  f_vcs->SetNpx(opts.n);
  fl->Close();
  delete fl;
  /**
   * Saving the convolution function to the output file
   */
  auto fl_out = TFile::Open(opts.ofname.c_str(), "recreate");
  fl_out->cd();
  f_vcs->Write();
  fl_out->Close();
  delete f_vcs;
  delete fl_out;
  return 0;
}

