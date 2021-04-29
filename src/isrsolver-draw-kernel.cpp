#include <iostream>
#include <string>
#include <functional>
#include <boost/program_options.hpp>
#include <TF1.h>
#include <TFile.h>
#include "KuraevFadin.hpp"
namespace po = boost::program_options;

/**
 * A part of program options
 */
typedef struct {
  /**
   * Number of points used to plot the function
   */
  int n;
  /**
   * The value of the center-of-mass energy
   */
  double energy;
  /**
   * Minimum x
   */
  double xmin;
  /**
   * Maximum x
   */
  double xmax;
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
      ("help,h", "message help")
      ("num-of-points,n", po::value<int>(&(opts->n))->default_value(10000), "number of center-of-mass energy points")
      ("energy,e", po::value<double>(&(opts->energy))->default_value(1.), "center-of-mass energy")
      ("xmin,m", po::value<double>(&(opts->xmin))->default_value(1.e-6), "minimum value of x")
      ("xmax,x", po::value<double>(&(opts->xmax))->default_value(0.5), "maximum value of x")
      ("ofname,o", po::value<std::string>(&(opts->ofname))->
       default_value("output_kuraev-fadin_kernel.root"),
       "output file path");
}

/**
 * Help message
 */
void help(const po::options_description& desc) {
  std::cout << desc << std::endl;
}

int main(int argc, char* argv[]) {
  po::options_description desc("   This tool is for drawing the kernel function F(x, s). Allowed options");
  CmdOptions opts;
  setOptions(&desc, &opts);
  po::variables_map vmap;
  po::store(po::parse_command_line(argc, argv, desc), vmap);
  po::notify(vmap);
  if (vmap.count("help")) {
    help(desc);
    return 0;
  }
  std::function<double(double*, double*)> ker_fcn =
      [opts](double* px, double*) {
        double result = kernelKuraevFadin(px[0], opts.energy * opts.energy);
        return result;
      };
  TF1 ker_f("kuraev_fadin_fcn", ker_fcn, opts.xmin, opts.xmax, 0);
  ker_f.SetNpx(opts.n);
  auto fl = TFile::Open(opts.ofname.c_str(), "recreate");
  fl->cd();
  ker_f.Write();
  fl->Close();
  delete fl;
  return 0;
}
