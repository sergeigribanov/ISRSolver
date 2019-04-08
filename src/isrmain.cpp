
#include <cstdlib>
#include <iostream>
#include <string>
#include <boost/program_options.hpp>
#include "TFile.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "RadSolver.h"

namespace po = boost::program_options;

struct cmdOptions {
  double thsd;
  std::string gname;
  std::string ifname;
  std::string ofname;
};

void setOptions(po::options_description* desc,
                struct cmdOptions* opts) {
  desc->add_options()
      ("help", "A simple tool, designed to find numerical"
       "solution of the Kuraev-Fadin equation.")
      ("thsd", po::value<double>(&(opts->thsd)), "Threshold (GeV).")
      ("gname", po::value<std::string>(&(opts->gname)),
       "Name of the measured cross section graph.")
      ("ifname", po::value<std::string>(&(opts->ifname)),
       "Path to input file.")
      ("ofname", po::value<std::string>(&(opts->ofname)),
       "Path to output file.");
}

void help(const po::options_description& desc) {
  std::cout << desc << std::endl;
}

int main(int argc, char* argv[]) {
  po::options_description desc("Allowed options:");
  struct cmdOptions opts;
  setOptions(&desc, &opts);
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  if (vm.count("help")) {
    help(desc);
    return 0;
  }
  if (vm.count("thsd") &&
      vm.count("gname") &&
      vm.count("ifname") &&
      vm.count("ofname")) {
    TMatrixT<double> invErrM;
    TMatrixT<double> integralOperatorMatrix;
    auto fl0 = TFile::Open(opts.ifname.c_str(), "read");
    auto gr = reinterpret_cast<TGraphErrors*>(fl0->Get(opts.gname.c_str()));
    RadSolver rs(gr, opts.thsd * opts.thsd);
    fl0->Close();
    auto gborn = rs.getBornCS(invErrM, integralOperatorMatrix);
    auto fl1 = TFile::Open(opts.ofname.c_str(), "recreate");
    fl1->cd();
    rs.visible_cs->Write("visible");
    gborn->Write("cborn");
    invErrM.Write("invBornErrMatrix");
    integralOperatorMatrix.Write("integralOperatorMatrix");
    fl1->Close();
  } else {
    help(desc);
  }
  return 0;
}
