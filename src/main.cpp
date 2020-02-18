#include <TAxis.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>

#include <boost/program_options.hpp>
#include <cstdlib>
#include <iostream>
#include <string>
#include "utils.hpp"
#include "RadSolver.hpp"
namespace po = boost::program_options;

typedef struct {
  double thsd;
  double spoint;
  std::string gname;
  std::string lbcs;
  std::string ifname;
  std::string ofname;
} CmdOptions;

void setOptions(po::options_description* desc, CmdOptions* opts) {
  desc->add_options()("help",
                      "A simple tool, designed to find numerical"
                      "solution of the Kuraev-Fadin equation.")(
      "thsd", po::value<double>(&(opts->thsd)), "Threshold (GeV).")(
      "gname", po::value<std::string>(&(opts->gname))->default_value("mcs"),
      "Name of the measured cross section graph.")(
      "lbcs", po::value<std::string>(&(opts->lbcs)),
      "Name of a TF1 object, which represents a left part of Born cross "
      "section "
      "if needed. The meaning of this function is the Bron cross section "
      "measured in "
      "previous experiments. Use this option, when you have Born cross section "
      "measurement, "
      "which is too far from threshold.")(
      "spoint", po::value<double>(&(opts->spoint)),
      "A value of a center-of-mass energy point, which is a start point for "
      "interpolation. "
      "Use this paramer in case if the left part of cross section is "
      "represented by the data "
      "obtained from the previous experiments.")(
      "ifname",
      po::value<std::string>(&(opts->ifname))->default_value("mcs.root"),
      "Path to input file.")(
      "ofname",
      po::value<std::string>(&(opts->ofname))->default_value("bcs.root"),
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
  auto measured_cs = find_object<TGraphErrors>(fl, opts.gname);
  RadSolver solver;
  solver.setThresholdEnergy(opts.thsd);
  solver.setMeasuredCrossSection(measured_cs);
  if (vmap.count("lbcs")) {
    if (vmap.count("spoint")) {
      solver.setStartPointEnergy(opts.spoint);
    }
    auto left_side_bcs = find_object<TF1>(fl, opts.lbcs);
    solver.setLeftSideOfBornCrossSection(left_side_bcs);
  }
  fl->Close();
  delete fl;
  solver.solve();
  solver.save(opts.ofname);
  return 0;
}
