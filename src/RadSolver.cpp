#include "RadSolver.h"
#include "sigmaCoefficients.h"
#include <algorithm>
#include "TVectorT.h"
#include <math.h>
#include <iostream>

RadSolver::RadSolver(TGraphErrors* gr, double s_treshold) {
  int N = gr->GetN();
  rightPart_.resize(N + 1);
  rightPart_[0].s_ = s_treshold;
  rightPart_[0].y_ = 0;
  rightPart_[0].ex_ = 0;
  rightPart_[0].ey_ = 0;
  for (int i = 0; i < N; ++i) {
    double ecm = gr->GetX() [i];
    rightPart_[i + 1].s_ = ecm * ecm;
    rightPart_[i + 1].y_ = gr->GetY()[i];
    rightPart_[i + 1].ex_ = gr->GetEX()[i];
    rightPart_[i + 1].ey_ = gr->GetEY()[i];
  }
  std::sort(rightPart_.begin(),  rightPart_.end(), [](const RightPart& x,
						      const RightPart& y) {
	      return x.s_ < y.s_;
	    });
  double X[N];
  double Y[N];
  double EX[N];
  double EY[N];
  
  for (int i = 0; i < N; ++i) {
    X[i] = sqrt(rightPart_[i+1].s_);
    Y[i] = rightPart_[i+1].y_;
    EX[i] = rightPart_[i+1].ex_;
    EY[i] = rightPart_[i+1].ey_;
  }
  visible_cs = new TGraphErrors (N, X, Y, 0, EY);
}

RadSolver::~RadSolver() {
  if (visible_cs != 0) {
    delete visible_cs;
  }
}

double RadSolver::getX(int n, int i) const {
  return 1 - rightPart_[i].s_ / rightPart_[n].s_;
}

TMatrixT<double> RadSolver::getEqMatrix () const {
  int N = rightPart_.size() - 1;
  TMatrixT<double> A (N, N);
  for (int n = 1; n <= N; ++n) {
    double sn = rightPart_[n].s_;
    for (int i = 1; i <= n; ++i) {
	double xm = getX (n, i - 1);
	double xi = getX(n, i);
	auto linc = getLinearSigmaCoeffs(xm, xi,
					 sn, xm, xi);
	A (n - 1 , i - 1) += linc [1];
	if (i > 1) {
	  A (n - 1 , i - 2) += linc [0];
	}	
    }
  }
  return (-1.)* A;
}

TGraphErrors* RadSolver::getBornCS(TMatrixT<double>& invErrM,
				   TMatrixT<double>& integralOperatorMatrix) const {
  auto eqM = getEqMatrix ();
  int N = rightPart_.size() - 1;
  TVectorT<double> ecm (N);
  TVectorT<double> ecm_err (N);
  TVectorT<double> cs (N);
  TVectorT<double> cs_err (N);
  TVectorT<double> vcs (N);
  TVectorT<double> vcs_err (N);

  for (int i = 0; i < N; ++i) {
    ecm (i) = sqrt (rightPart_[i + 1].s_);
    ecm_err (i) = rightPart_[i + 1].ex_;
    vcs (i) =  rightPart_[i + 1].y_;
    vcs_err (i) = rightPart_[i + 1].ey_;
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

  invErrM.ResizeTo (N, N);
  invErrM = eqMT * Lam * eqM;

  integralOperatorMatrix.ResizeTo (N, N);
  integralOperatorMatrix = eqM;

  
  TMatrixT<double> errM (N, N);

  errM = invErrM;

  errM.Invert();

  for (int i = 0; i < N; ++i) {
    cs_err (i) = sqrt(errM (i, i));
  }

  TGraphErrors* gr = new TGraphErrors(N, ecm.GetMatrixArray(),
                                      cs.GetMatrixArray(),
                                      ecm_err.GetMatrixArray(),
                                      cs_err.GetMatrixArray());
  
  return gr;
}
