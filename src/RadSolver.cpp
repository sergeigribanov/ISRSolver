#include <math.h>
#include <iostream>
#include <algorithm>
#include <TVectorT.h>
#include <TFile.h>
#include "RadSolver.h"
#include "sigmaCoefficients.h"

RadSolver::RadSolver() : _s_threshold(0) {
}

RadSolver::~RadSolver() {
}

double RadSolver::getX(int n, int i) const {
  return 1 - _measured_cs_data[i].s_ / _measured_cs_data[n].s_;
}

TMatrixT<double> RadSolver::getEqMatrix () const {
  int N = _measured_cs_data.size() - 1;
  TMatrixT<double> A (N, N);
  for (int n = 1; n <= N; ++n) {
    double sn = _measured_cs_data[n].s_;
    for (int i = 1; i <= n; ++i) {
	double xm = getX (n, i - 1);
	double xi = getX(n, i);
	auto linc = getLinearSigmaCoeffs(xm, xi, sn, xm, xi);
	A (n - 1 , i - 1) += linc [1];
	if (i > 1) {
	  A (n - 1 , i - 2) += linc [0];
	}	
    }
  }
  return (-1.)* A;
}

void RadSolver::solve() {
  auto eqM = getEqMatrix ();
  int N = _measured_cs_data.size() - 1;
  TVectorT<double> ecm (N);
  TVectorT<double> ecm_err (N);
  TVectorT<double> cs (N);
  TVectorT<double> cs_err (N);
  TVectorT<double> vcs (N);
  TVectorT<double> vcs_err (N);

  for (int i = 0; i < N; ++i) {
    ecm (i) = sqrt (_measured_cs_data[i + 1].s_);
    ecm_err (i) = _measured_cs_data[i + 1].ex_;
    vcs (i) =  _measured_cs_data[i + 1].y_;
    vcs_err (i) = _measured_cs_data[i + 1].ey_;
  }

  TMatrixT<double> eqMT (N, N);
  eqMT.Transpose (eqM);
  
  TMatrixT<double> eqMI (N, N);
  eqMI = eqM;
  eqMI.Invert ();

  cs = eqMI * vcs;

  TMatrixT<double> Lam (N, N);
  for (int i = 0; i < N; ++i) {
    Lam (i, i) = 1. / vcs_err (i) / vcs_err (i);
  }

  _inverse_error_matrix.ResizeTo (N, N);
  _inverse_error_matrix = eqMT * Lam * eqM;

  _integral_operator_matrix.ResizeTo (N, N);
  _integral_operator_matrix = eqM;

  
  TMatrixT<double> errM (N, N);

  errM = _inverse_error_matrix;

  errM.Invert();

  for (int i = 0; i < N; ++i) {
    cs_err (i) = sqrt(errM (i, i));
  }

  _born_cs =  TGraphErrors(N, ecm.GetMatrixArray(),
			   cs.GetMatrixArray(),
			   ecm_err.GetMatrixArray(),
			   cs_err.GetMatrixArray());
}

void RadSolver::setMeasuredCrossSection(TGraphErrors* graph) {
  const int N = graph->GetN();
  _measured_cs_data.resize(N + 1);
  _measured_cs_data[0].s_ = _s_threshold;
  _measured_cs_data[0].y_ = 0;
  _measured_cs_data[0].ex_ = 0;
  _measured_cs_data[0].ey_ = 0;
  for (int i = 0; i < N; ++i) {
    double ecm = graph->GetX() [i];
    _measured_cs_data[i + 1].s_ = ecm * ecm;
    _measured_cs_data[i + 1].y_ = graph->GetY()[i];
    _measured_cs_data[i + 1].ex_ = graph->GetEX()[i];
    _measured_cs_data[i + 1].ey_ = graph->GetEY()[i];
  }
  std::sort(_measured_cs_data.begin(),
	    _measured_cs_data.end(),
	    [](const RightPart& x, const RightPart& y)
	    {return x.s_ < y.s_;});
  _measured_cs.Set(0);
  for (int i = 0; i < N; ++i) {
    _measured_cs.SetPoint(i, sqrt(_measured_cs_data[i+1].s_),
			  _measured_cs_data[i+1].y_);
    _measured_cs.SetPointError(i, _measured_cs_data[i+1].ex_,
			       _measured_cs_data[i+1].ey_);
  }
}

const TGraphErrors& RadSolver::getBornCrossSection() const {
  return _born_cs;
}

const TGraphErrors& RadSolver::getMeasuredCrossSection() const {
  return _measured_cs;
}

const TMatrixT<double>& RadSolver::getIntegralOeratorMatrix() const {
  return _integral_operator_matrix;
}

const TMatrixT<double>& RadSolver::getInverseErrorMatrix() const {
  return _inverse_error_matrix;
}

void RadSolver::save(const std::string& path) {
  auto fl = TFile::Open(path.c_str(), "recreate");
  fl->cd();
  _measured_cs.Write("measured_cs");
  _born_cs.Write("born_cs");
  _integral_operator_matrix.Write("integral_operator_matrix");
  _inverse_error_matrix.Write("inverse_error_matrix");
  fl->Close();
  delete fl;
}

double RadSolver::getThresholdS() const {
  return _s_threshold;
}

void RadSolver::setThresholdS(double s_threshold) {
  _s_threshold = s_threshold;
}
