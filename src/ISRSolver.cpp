#include "ISRSolver.hpp"

#include <TFile.h>
#include <TF1.h>
#include <TGraph.h>

#include <functional>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <Eigen/Core>

#include "kuraev_fadin.hpp"

ISRSolver::ISRSolver(const std::string& inputPath,
		     const InputOptions& inputOpts) :
  _inputOpts(inputOpts),
  _sT(inputOpts.thresholdEnergy * inputOpts.thresholdEnergy) {
  auto fl = TFile::Open(inputPath.c_str(), "read");
  auto graph = dynamic_cast<TGraphErrors*>
    (fl->Get(_inputOpts.measuredCSGraphName.c_str()));
  _n = graph->GetN();
  double energyKoeff; 
  if (_inputOpts.energyUnitMeVs)
    energyKoeff = 1.e-3;
  else
    energyKoeff = 1;
  std::vector<CSData> measuredCS;
  measuredCS.reserve(_n);
  for (std::size_t i = 0; i < _n; ++i)
    measuredCS.push_back
      ({
	.cmEnergy = graph->GetX()[i] * energyKoeff,
	.cs = graph->GetY()[i],
	.cmEnergyError = graph->GetEX()[i] * energyKoeff,
	.csError = graph->GetEY()[i]});
  fl->Close();
  std::sort(measuredCS.begin(), measuredCS.end(),
	    [](const CSData& x, const CSData& y)
	    {return x.cmEnergy < y.cmEnergy;});
  _measuredCSData = {
    .cmEnergy = Eigen::VectorXd(_n),
    .s = Eigen::VectorXd(_n),
    .cs = Eigen::VectorXd(_n),
    .cmEnergyError = Eigen::VectorXd(_n),
    .invCSErrMatrix = Eigen::MatrixXd::Zero(_n, _n)};
  std::transform(measuredCS.begin(),
		 measuredCS.end(),
		 _measuredCSData.cmEnergy.data(),
		 [](const CSData& x)
		 {return x.cmEnergy;});
  std::transform(measuredCS.begin(),
		 measuredCS.end(),
		 _measuredCSData.s.data(),
		 [](const CSData& x)
		 {return x.cmEnergy * x.cmEnergy;});
  std::transform(measuredCS.begin(),
		 measuredCS.end(),
		 _measuredCSData.cs.data(),
		 [](const CSData& x)
		 {return x.cs;});
  std::transform(measuredCS.begin(),
		 measuredCS.end(),
		 _measuredCSData.cmEnergyError.data(),
		 [](const CSData& x)
		 {return x.cmEnergyError;});
  Eigen::VectorXd tmpv(_n);
  std::transform(measuredCS.begin(),
		 measuredCS.end(),
		 tmpv.data(),
		 [](const CSData& x) {return 1. / x.csError / x.csError;});
  _measuredCSData.invCSErrMatrix.diagonal() = tmpv;

  _xPoints.resize(_n + 1);
  _xPoints[0] = 0;
  double* sVec = _measuredCSData.s.data();
  std::transform(sVec, sVec + getN(),
		 _xPoints.begin() + 1,
		 [this](double s)
		 {return 1 - this->_sT / s;});
  
  setDefaultInterpSettings();
}

ISRSolver::~ISRSolver() {}

void ISRSolver::testPrint() const {
  Eigen::VectorXd y(getN());
  for (std::size_t i = 0; i < _n; ++i) {
    y(i) = std::sin(10. * (_measuredCSData.cmEnergy(i) - _inputOpts.thresholdEnergy));
  }
  std::function<double(double*, double*)> fcn = 
    [this, y](double* x, double* par) {
      double en = x[0];
      if (en <= _inputOpts.thresholdEnergy)
	return 0.;
      std::size_t i = 0;
      double enp = _inputOpts.thresholdEnergy;
      while (i < this->getN()) { 
	if (en > enp &&
	    en <= _measuredCSData.cmEnergy(i))
	  break;
	enp = _measuredCSData.cmEnergy(i);
	i++;
      }
      if (i == this->getN())
	i--;
      
      Eigen::VectorXd coeffs =
	this->interpInvMatrix(i) *
	this->permutation(i) * y;
      double result = 0;
      for (int k = 0; k < coeffs.size(); ++k)
	result += coeffs(k) * std::pow(en * en, k);
      return result;
    };
  auto f1 = new TF1("f1", fcn, 1.19, 2., 0);
  f1->SetNpx(1000);
  TGraph gr(this->getN(), _measuredCSData.cmEnergy.data(), y.data());
  auto fl = TFile::Open("test.root", "recreate");
  fl->cd();
  f1->Write("f1");
  gr.Write("gr");
  fl->Close();
}

double ISRSolver::getXmin(int j) const {
  return _xPoints[j];
}

double ISRSolver::getXmax(int j) const {
  return _xPoints[j + 1];
}

std::size_t ISRSolver::getN() const {
  return _n;
}

void ISRSolver::setDefaultInterpSettings() {
  _interpSettings.resize(_n, {.numCoeffs = 5, .nPointsLeft = 2});
  _interpSettings[0] = {.numCoeffs = 2, .nPointsLeft = 1};
  _interpSettings[_n - 2] = {.numCoeffs = 3, .nPointsLeft = 1};
  _interpSettings[_n - 1] = {.numCoeffs = 2, .nPointsLeft = 1};
}

double ISRSolver::getNumCoeffs(int j) const {
  return _interpSettings[j].numCoeffs;
}

double ISRSolver::getNPointsLeft(int j) const {
  return _interpSettings[j].nPointsLeft;
}

Eigen::MatrixXd ISRSolver::permutation(int j) const {
  const std::size_t nc = getNumCoeffs(j);
  const std::size_t npl = getNPointsLeft(j);
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(nc, getN());
  int p = j - npl;
  for (std::size_t l = 0; l < nc; ++l) {
    if (p >= 0) {
      result(l, p) = 1;
    }
    p++;
  }
  return result;
}

Eigen::MatrixXd ISRSolver::interpInvMatrix(int j) const {
  const std::size_t nc = getNumCoeffs(j);
  const std::size_t npl = getNPointsLeft(j);
  Eigen::MatrixXd interpMatrix = Eigen::MatrixXd::Zero(nc, nc);
  int p = j - npl;
  for (std::size_t l = 0; l < nc; ++l) {
    if (p < 0) {
      for (std::size_t k = 0; k < nc; ++k)
	interpMatrix(l, k) = std::pow(_sT, k);
    } else {
      for (std::size_t k = 0; k < nc; ++k)
	interpMatrix(l, k) = std::pow(_measuredCSData.s(p), k);
    }
     p++;
  }
  return interpMatrix.inverse();
}

Eigen::MatrixXd ISRSolver::polConvKuraevFadinMatrix(int j) const {
  const std::size_t nc = getNumCoeffs(j);
  Eigen::MatrixXd result(getN(), nc);
  std::size_t k;
  std::function<double(double)> fcn =
    [&k](double x) {
      return std::pow(1 - x, k);
    };
  for (std::size_t i = 0; i < getN(); ++i)
    for (k = 0; k < nc; ++k) {
      result(i, k) = std::pow(_measuredCSData.s(i), k) *
	kuraev_fadin_convolution(_measuredCSData.s(i), fcn, getXmin(j), getXmax(j));
    }
  return result;
}

Eigen::MatrixXd ISRSolver::evalA(int j) const {
  return polConvKuraevFadinMatrix(j) * interpInvMatrix(j) * permutation(j);
}

Eigen::MatrixXd ISRSolver::evalEqMatrix() const {
  // !!! TO DO
}

// void ISRSolver::solve() { 
//   auto eqM = getEqMatrix();
//   int N = _measured_cs_data.size() - 1;
//   Eigen::VectorXd ecm = Eigen::VectorXd::Zero(N);
//   Eigen::VectorXd ecm_err = Eigen::VectorXd::Zero(N);
//   Eigen::VectorXd mcs = Eigen::VectorXd::Zero(N);
//   Eigen::VectorXd mcs_err = Eigen::VectorXd::Zero(N);
//   Eigen::VectorXd cs;
//   Eigen::VectorXd cs_err;
//   std::function<double(double)> lsbcs =
//     [this](double s)
//     {
//       double e0 = this->getThresholdEnergy();
//       double e1 = this->_left_side_bcs->GetXmin();
//       double en = std::sqrt(s);
//       if (e0 < e1 && en < e1) { 
// 	return this->_left_side_bcs->Eval(e1) * (en - e0) / (e1 - e0);
//       }
//       return this->_left_side_bcs->Eval(en);
//   };
//   double s_start = _start_point_energy * _start_point_energy;
//   double s_threshold = _threshold_energy * _threshold_energy;
//   for (int i = 0; i < N; ++i) {
//     ecm(i) = sqrt(_measured_cs_data[i + 1].s);
//     ecm_err(i) = _measured_cs_data[i + 1].ex;
//     mcs(i) = _measured_cs_data[i + 1].y;
//     if (_left_side_bcs) {
//       mcs(i) -= kuraev_fadin_convolution(
//           _measured_cs_data[i + 1].s, lsbcs,
// 	  1 - s_start / _measured_cs_data[i + 1].s,
//           1 - s_threshold / _measured_cs_data[i + 1].s);
//     }
//     mcs_err(i) = _measured_cs_data[i + 1].ey;
//   }
//   cs = eqM.completeOrthogonalDecomposition().solve(mcs);
//   Eigen::MatrixXd lam = Eigen::MatrixXd::Zero(N, N);
//   lam.diagonal() = mcs_err.array().square().inverse();
//   Eigen::MatrixXd invErrM;
//   invErrM = eqM.transpose() * lam * eqM;
//   _inverse_error_matrix.ResizeTo(N, N);
//   _integral_operator_matrix.ResizeTo(N, N);
//   _integral_operator_matrix.SetMatrixArray(eqM.data());
//   _integral_operator_matrix.Transpose(_integral_operator_matrix);
//   _inverse_error_matrix.SetMatrixArray(invErrM.data());
//   _inverse_error_matrix.Transpose(_inverse_error_matrix);
//   cs_err = invErrM.inverse().diagonal().array().sqrt();
//   _born_cs =
//       TGraphErrors(N, ecm.data(), cs.data(), ecm_err.data(), cs_err.data());
// }


// void ISRSolver::save(const std::string& path) {
//   auto fl = TFile::Open(path.c_str(), "recreate");
//   fl->cd();
//   _measured_cs.Write("measured_cs");
//   _born_cs.Write("born_cs");
//   _integral_operator_matrix.Write("integral_operator_matrix");
//   _inverse_error_matrix.Write("inverse_error_matrix");
//   fl->Close();
//   delete fl;
// }



// std::pair<double, double> ISRSolver::coeffs(double xm, double xi, double s) {
//   double det = xm - xi;
//   double integral0 = kuraev_fadin_polinomial_convolution(s, xm, xi, 0);
//   double integral1 = kuraev_fadin_polinomial_convolution(s, xm, xi, 1);
//   double cm = (integral1 - xi * integral0) / det;
//   double ci = (-integral1 + xm * integral0) / det;
//   return std::make_pair(cm, ci);
// }
