#include <iostream>
#include <algorithm>
#include <functional>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include "Interpolator.hpp"

void test_basis(const Eigen::VectorXd& y,
                const Interpolator& interp,
                double energyMin = 0.827,
                double energyMax = 2.5) {
  std::function<double(double*, double*)> fcn =
      [interp, y](double* px, double*) {
        double result = 0;
        for (int i = 0; i < y.rows(); ++i) {
          result += interp.basisEval(i, px[0]) * y(i);
        }
        return result;
      };
  std::function<double(double*, double*)> dfcn =
      [interp, y](double* px, double*) {
        double result = 0;
        for (int i = 0; i < y.rows(); ++i) {
          result += interp.basisDerivEval(i, px[0]) * y(i);
        }
        return result;
      };
  TF1 f0("f_basis", fcn, energyMin, energyMax, 0);
  f0.SetNpx(10000);
  TF1 f0_deriv("f_deriv_basis", dfcn, energyMin, energyMax, 0);
  f0_deriv.SetNpx(10000);
  TFile fl("basis_test.root", "recreate");
  fl.cd();
  f0.Write();
  f0_deriv.Write();
  fl.Close();
}

void test_entire(const Eigen::VectorXd& y,
                 const Interpolator& interp,
                 double energyMin = 0.827,
                 double energyMax = 2.5) {
  std::function<double(double*, double*)> fcn =
      [interp, y](double* px, double*) {
        return interp.eval(y, px[0]);
      };
  std::function<double(double*, double*)> dfcn =
      [interp, y](double* px, double*) {
        return interp.derivEval(y, px[0]);
      };
  TF1 f0("f_entire", fcn, energyMin, energyMax, 0);
  f0.SetNpx(10000);
  TF1 f0_deriv("f_deriv_entire", dfcn, energyMin, energyMax, 0);
  f0_deriv.SetNpx(10000);
  TFile fl("entire_test.root", "recreate");
  fl.cd();
  f0.Write();
  f0_deriv.Write();
  fl.Close();
}

int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::cout << "[!] Wrong number of command-line arguments" << std::endl;
  }
  auto fl = TFile::Open(argv[1], "read");
  auto bcs = dynamic_cast<TGraphErrors*>(fl->Get("bcs"));
  const int n = bcs->GetN();
  std::vector<double> x(n);
  std::vector<double> y(n);
  std::copy(bcs->GetX(), bcs->GetX() + n, x.begin());
  std::copy(bcs->GetY(), bcs->GetY() + n, y.begin());
  fl->Close();
  const double thresholdEnergy = 0.827;
  Eigen::Map<Eigen::VectorXd> cmEnergies(x.data(), x.size());

  std::vector<std::tuple<bool, int, int>> interpRangeSettings;
  interpRangeSettings.push_back(std::make_tuple(false, 0, 5));
  interpRangeSettings.push_back(std::make_tuple(true, 6, 40));
  interpRangeSettings.push_back(std::make_tuple(false, 41, 49));
  Interpolator interp(interpRangeSettings, cmEnergies, thresholdEnergy);
 
  // Interpolator interp(cmEnergies, thresholdEnergy);
  Eigen::Map<Eigen::VectorXd> v(y.data(), y.size());
  test_basis(v, interp);
  test_entire(v, interp);
  return 0;
}
