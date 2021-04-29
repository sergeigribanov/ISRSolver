#include <iostream>
#include <vector>
#include <functional>
#include <boost/program_options.hpp>
#include <TFile.h>
#include <TF1.h>
#include <TEfficiency.h>
#include "KuraevFadin.hpp"
namespace po = boost::program_options;

/**
 * A part of program options
 */
typedef struct {
  /**
   * Number of points used to plot radiative
   * correction function
   */
  std::size_t n;
  /**
   * Threshold energy
   */
  double thsd;
  /**
   * Minimum center-of-mass energy
   */
  double minen;
  /**
   * Maximum center-of-mass energy
   */
  double maxen;
  /**
   * Path to the input .root file that contains
   * Born cross section function
   */
  std::string ifname;
  /**
   * Name of the Born cross section function
   */
  std::string bcs_fcn_name;
  /**
   * Name of the detection efficiency object (TEfficiency)
   */
  std::string efficiency_name;
  /**
   * Path to the output file
   */
  std::string ofname;
} CmdOptions;

/**
 * Setting up program options
 */
void setOptions(po::options_description* desc, CmdOptions* opts) {
  desc->add_options()
      ("help,h", "radiative correction calculator")
      ("thsd,t", po::value<double>(&(opts->thsd)), "threshold (GeV)")
      ("number-of-radocrr-points,n", po::value<std::size_t>(&(opts->n))->default_value(1000),
       "number of points to tabulate ratiative correction (delta)")
      ("minen,m",  po::value<double>(&(opts->minen)), "minimum energy")
      ("maxen,x",  po::value<double>(&(opts->maxen)), "maximum energy")

      ( "ifname,i",
        po::value<std::string>(&(opts->ifname))->default_value("input.root"),
        "path to input file")
      ("ofname,o",
       po::value<std::string>(&(opts->ofname))->default_value("output.root"),
       "path to output file")
      ("born-cs-fcn-name,b",
       po::value<std::string>(&(opts->bcs_fcn_name))->default_value("f_bcs"),
       "the name of the Born cross section function (TF1*)")
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
  TEfficiency* teff = nullptr;
  /**
   * Opening input file
   */
  auto fl = TFile::Open(opts.ifname.c_str(), "read");
  /**
   * Reading the Born cross section function from the input file
   */
  auto fbcs = dynamic_cast<TF1*>(fl->Get(opts.bcs_fcn_name.c_str())->Clone());
  /**
   * Reading the detection efficiency from the input file
   */
  if (vmap.count("efficiency-name")) {
    teff = dynamic_cast<TEfficiency*>(fl->Get(opts.efficiency_name.c_str())->Clone());
  }
  fl->Close();
  delete fl;
  /**
   * Converting the Born cross section function to
   * a form of std::function
   */
  std::function<double(double)> fcn =
      [fbcs](double energy) {
        double result = fbcs->Eval(energy);
        return result;
      };
  /**
   * Converting the detection efficiency to a form
   * of std::function
   */
  std::function<double(double, double)> eff =
      [teff](double x, double en) {
        int bin = teff->FindFixBin(x, en);
        double result = teff->GetEfficiency(bin);
        return result;
      };
  std::vector<double> ens;
  std::vector<double> radcorrs;
  ens.reserve(opts.n);
  radcorrs.reserve(opts.n);
  const double s_th = opts.thsd * opts.thsd;
  /**
   * Creating radiative correction function
   */
  std::function<double(double*, double*)> radcorr_fcn =
      [s_th, teff,fcn, eff](double* px, double*) {
        const double en = px[0];
        const double fe = fcn(en);
        if (fe < 1.e-12) {
          return 0.;
        }
        if (teff) {
          double result = convolutionKuraevFadin(en, fcn, 0, 1 - s_th / en / en, eff) / fe - 1;
          return result;
        }
        double result = convolutionKuraevFadin(en, fcn, 0, 1 - s_th / en / en) / fe - 1;
        return result;
      };
  TF1 radf("radcorr", &radcorr_fcn, opts.minen, opts.maxen, 0);
  radf.SetNpx(opts.n);
  /**
   * Saving radiative correction function
   * to the output file
   */
  auto ofl = TFile::Open(opts.ofname.c_str(), "recreate");
  ofl->cd();
  radf.Write();
  delete ofl;
  delete fbcs;
  delete teff;
  return 0;
}
