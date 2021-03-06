#include <iostream>
#include <algorithm>
#include <string>
#include <boost/program_options.hpp>
#include <TFile.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include "Interpolator.hpp"
namespace po = boost::program_options;

/**
 * A part of program options
 */
typedef struct {
  /**
   * Threshold energy
   */
  double thsd;
  /**
   * Name of the graph (TGraphErrors) to be interpolated
   */
  std::string graph_name;
  /**
   * Input file path
   */
  std::string ifname;
  /**
   * Output file path
   */
  std::string ofname;
  /**
   * Path to the file that contains interpolation
   * settings
   */
  std::string interp;
} CmdOptions;

/**
 * Setting up program options
 */
void setOptions(po::options_description* desc, CmdOptions* opts) {
  desc->add_options()
      ("help,h", "help message")
      ("thsd,t", po::value<double>(&(opts->thsd)), "threshold (GeV)")
      ("cs-graph-name,c",
       po::value<std::string>(&(opts->graph_name))->default_value("bcs"),
       "name of a cross section graph (TGraphErrors*)")
      ("ifname,i",
       po::value<std::string>(&(opts->ifname))->default_value("input.root"),
       "path to input file")
      ("ofname,o",
       po::value<std::string>(&(opts->ofname))->default_value("output.root"),
       "path to output file")
      ("interp,r",
       po::value<std::string>(&(opts->interp)),
       "path to JSON file with interpolation settings");
}

/**
 * Help message
 */
void help(const po::options_description& desc) {
  std::cout << desc << std::endl;
}

/**
 * Interpolation
 */
void interpBasis(const Eigen::VectorXd& y,
                 const Interpolator& interp,
                 double emax,
                 const CmdOptions& opts) {
  /**
   * Creating interpolation function
   */
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
  /**
   * Saving the interpolation function to the output
   * file
   */
  auto fl = TFile::Open(opts.ofname.c_str(), "recreate");
  fl->cd();
  f0.Write();
  fl->Close();
  delete fl;
}

int main(int argc, char* argv[]) {
  po::options_description desc("   This tool is for testing cross section interpolation. Allowed options");
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
