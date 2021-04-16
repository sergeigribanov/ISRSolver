#include <iostream>
#include <string>
#include <functional>
#include <boost/program_options.hpp>
#include <TF1.h>
#include <TFile.h>
#include "kuraev_fadin.hpp"
namespace po = boost::program_options;

typedef struct {
  int n;
  double energy;
  double xmin;
  double xmax;
  std::string ofname;
} CmdOptions;

void setOptions(po::options_description* desc, CmdOptions* opts) {
  desc->add_options()
      ("help,h", "Drawing Integral of Kuraev-Fadin kernel function")
      ("num-of-points,n", po::value<int>(&(opts->n))->default_value(10000), "Center-of-mass energy")
      ("energy,e", po::value<double>(&(opts->energy))->default_value(1.), "Center-of-mass energy")
      ("xmin,m", po::value<double>(&(opts->xmin))->default_value(1.e-6), "Minimum value of x")
      ("xmax,x", po::value<double>(&(opts->xmax))->default_value(0.5), "Maximum value of x")
      ("ofname,o", po::value<std::string>(&(opts->ofname))->
       default_value("output_kuraev-fadin_kernel_integral.root"),
       "Output file path");
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
  std::function<double(double*, double*)> ker_int_fcn =
      [opts](double* px, double*) {
        double result = kuraev_fadin_convolution(
            opts.energy, [](double) {return 1.;}, 0, px[0]);
        return result;
      };
  TF1 ker_int_f("kuraev_fadin_integral_fcn", ker_int_fcn, opts.xmin, opts.xmax, 0);
  ker_int_f.SetNpx(opts.n);
  auto fl = TFile::Open(opts.ofname.c_str(), "recreate");
  fl->cd();
  ker_int_f.Write();
  fl->Close();
  delete fl;
  return 0;
}
