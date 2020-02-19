#include "ISRSolver.hpp"

#include <TFile.h>

#include <algorithm>
#include <cmath>
#include <iostream>

#include "kuraev_fadin.hpp"

ISRSolver::ISRSolver()
    : _threshold(false),
      _start_point(false),
      _threshold_energy(0),
      _start_point_energy(0),
      _left_side_bcs(nullptr) {}

ISRSolver::~ISRSolver() {
  if (_left_side_bcs) {
    delete _left_side_bcs;
  }
}

double ISRSolver::getX(int n, int i) const {
  return 1 - _measured_cs_data[i].s / _measured_cs_data[n].s;
}

Eigen::MatrixXd ISRSolver::getEqMatrix() const {
  int N = _measured_cs_data.size() - 1;
  Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(N, N);
  for (int n = 1; n <= N; ++n) {
    double sn = _measured_cs_data[n].s;
    for (int i = 1; i <= n; ++i) {
      double xm = getX(n, i - 1);
      double xi = getX(n, i);
      auto linc = coeffs(xm, xi, sn);
      matrix(n - 1, i - 1) += linc.second;
      if (i > 1) {
        matrix(n - 1, i - 2) += linc.first;
      }
    }
  }
  return -matrix;
}

void ISRSolver::solve() {
  check();
  auto eqM = getEqMatrix();
  int N = _measured_cs_data.size() - 1;
  Eigen::VectorXd ecm = Eigen::VectorXd::Zero(N);
  Eigen::VectorXd ecm_err = Eigen::VectorXd::Zero(N);
  Eigen::VectorXd mcs = Eigen::VectorXd::Zero(N);
  Eigen::VectorXd mcs_err = Eigen::VectorXd::Zero(N);
  Eigen::VectorXd cs;
  Eigen::VectorXd cs_err;
  std::function<double(double)> lsbcs =
    [this](double s)
    {
      double e0 = this->getThresholdEnergy();
      double e1 = this->_left_side_bcs->GetXmin();
      double en = std::sqrt(s);
      if (e0 < e1 && en < e1) { 
	return this->_left_side_bcs->Eval(e1) * (en - e0) / (e1 - e0);
      }
      return this->_left_side_bcs->Eval(en);
  };
  double s_start = _start_point_energy * _start_point_energy;
  double s_threshold = _threshold_energy * _threshold_energy;
  for (int i = 0; i < N; ++i) {
    ecm(i) = sqrt(_measured_cs_data[i + 1].s);
    ecm_err(i) = _measured_cs_data[i + 1].ex;
    mcs(i) = _measured_cs_data[i + 1].y;
    if (_left_side_bcs) {
      mcs(i) -= kuraev_fadin_convolution(
          _measured_cs_data[i + 1].s, lsbcs,
	  1 - s_start / _measured_cs_data[i + 1].s,
          1 - s_threshold / _measured_cs_data[i + 1].s);
    }
    mcs_err(i) = _measured_cs_data[i + 1].ey;
  }
  cs = eqM.completeOrthogonalDecomposition().solve(mcs);
  Eigen::MatrixXd lam = Eigen::MatrixXd::Zero(N, N);
  lam.diagonal() = mcs_err.array().square().inverse();
  Eigen::MatrixXd invErrM;
  invErrM = eqM.transpose() * lam * eqM;
  _inverse_error_matrix.ResizeTo(N, N);
  _integral_operator_matrix.ResizeTo(N, N);
  _integral_operator_matrix.SetMatrixArray(eqM.data());
  _integral_operator_matrix.Transpose(_integral_operator_matrix);
  _inverse_error_matrix.SetMatrixArray(invErrM.data());
  _inverse_error_matrix.Transpose(_inverse_error_matrix);
  cs_err = invErrM.inverse().diagonal().array().sqrt();
  _born_cs =
      TGraphErrors(N, ecm.data(), cs.data(), ecm_err.data(), cs_err.data());
}

void ISRSolver::setMeasuredCrossSection(const TGraphErrors* graph) {
  const int N = graph->GetN();
  _measured_cs_data.resize(N + 1);
  for (int i = 0; i < N; ++i) {
    double ecm = graph->GetX()[i];
    _measured_cs_data[i + 1].s = ecm * ecm;
    _measured_cs_data[i + 1].y = graph->GetY()[i];
    _measured_cs_data[i + 1].ex = graph->GetEX()[i];
    _measured_cs_data[i + 1].ey = graph->GetEY()[i];
  }
  std::sort(_measured_cs_data.begin(), _measured_cs_data.end(),
            [](const RightPart& x, const RightPart& y) { return x.s < y.s; });
  _measured_cs.Set(0);
  for (int i = 0; i < N; ++i) {
    _measured_cs.SetPoint(i, sqrt(_measured_cs_data[i + 1].s),
                          _measured_cs_data[i + 1].y);
    _measured_cs.SetPointError(i, _measured_cs_data[i + 1].ex,
                               _measured_cs_data[i + 1].ey);
  }
}

void ISRSolver::setLeftSideOfBornCrossSection(const TF1* lbcs) {
  if (_left_side_bcs) {
    delete _left_side_bcs;
  }
  _left_side_bcs = dynamic_cast<TF1*>(lbcs->Clone("left_side_bcs"));
}

const TGraphErrors& ISRSolver::getBornCrossSection() const { return _born_cs; }

const TGraphErrors& ISRSolver::getMeasuredCrossSection() const {
  return _measured_cs;
}

const TMatrixT<double>& ISRSolver::getIntegralOeratorMatrix() const {
  return _integral_operator_matrix;
}

const TMatrixT<double>& ISRSolver::getInverseErrorMatrix() const {
  return _inverse_error_matrix;
}

void ISRSolver::save(const std::string& path) {
  auto fl = TFile::Open(path.c_str(), "recreate");
  fl->cd();
  _measured_cs.Write("measured_cs");
  _born_cs.Write("born_cs");
  _integral_operator_matrix.Write("integral_operator_matrix");
  _inverse_error_matrix.Write("inverse_error_matrix");
  fl->Close();
  delete fl;
}

double ISRSolver::getThresholdEnergy() const { return _threshold_energy; }

void ISRSolver::setThresholdEnergy(double energy) {
  _threshold_energy = energy;
  _threshold = true;
}

bool ISRSolver::isThresholdSEnabled() const { return _threshold; }

bool ISRSolver::isStartSEnabled() const { return _start_point; }

void ISRSolver::disableThreshold() { _threshold = false; }

void ISRSolver::enableThreshold() { _threshold = true; }

void ISRSolver::disableStartPoint() { _start_point = false; }

void ISRSolver::enableStartPoint() { _start_point = true; }

void ISRSolver::setStartPointEnergy(double energy) {
  _start_point_energy = energy;
  _start_point = true;
}

void ISRSolver::check() {
  if (_threshold && _measured_cs.GetX()[0] <= _threshold_energy) {
    std::cerr << "[!] Energies can not be equal or lower than threshold."
              << std::endl;
    exit(1);
  }
  if (_left_side_bcs) {
    if (!_threshold) {
      _threshold_energy = _left_side_bcs->GetXmin();
    } else {
      if (_threshold_energy >= _left_side_bcs->GetXmax()) {
        std::cerr
            << "[!] Threshold energy can not be equal or higher "
               "than maximum X value of left side Born cross section function."
            << std::endl;
        exit(1);
      }
    }
    if (!_start_point) {
      _start_point_energy = _left_side_bcs->GetXmax();
    }
    if (_start_point_energy <= _threshold_energy) {
      std::cerr << "[!] Start energy can not be equal or lower than threshold."
                << std::endl;
      exit(1);
    }
    if (_start_point_energy > _left_side_bcs->GetXmax()) {
      std::cerr
          << "[!] Start energy can not be higher than threshold."
          << "than maximum X value of left side Born cross section function."
          << std::endl;
      exit(1);
    }
  }
  if (!_threshold && !_left_side_bcs) {
    std::cerr << "[!] You need to set a threshold energy" << std::endl;
    exit(1);
  }
  if (_threshold && !_left_side_bcs) {
    _start_point_energy = _threshold_energy;
  }
  _measured_cs_data[0].s = _start_point_energy * _start_point_energy;
  if (_left_side_bcs) {
    _measured_cs_data[0].y = _left_side_bcs->Eval(_start_point_energy);
  } else {
    _measured_cs_data[0].y = 0;
  }
  _measured_cs_data[0].ex = 0;
  _measured_cs_data[0].ey = 0;
}

std::pair<double, double> ISRSolver::coeffs(double xm, double xi, double s) {
  double det = xm - xi;
  double integral0 = kuraev_fadin_polinomial_convolution(s, xm, xi, 0);
  double integral1 = kuraev_fadin_polinomial_convolution(s, xm, xi, 1);
  double cm = (integral1 - xi * integral0) / det;
  double ci = (-integral1 + xm * integral0) / det;
  return std::make_pair(cm, ci);
}
