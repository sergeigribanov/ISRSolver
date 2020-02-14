#ifndef __RADSOLVER_HPP__
#define __RADSOLVER_HPP__
#include <TF1.h>
#include <TGraphErrors.h>
#include <TMatrixT.h>

#include <Eigen/Dense>
#include <string>
#include <utility>
#include <vector>

struct RightPart {
  double s;
  double y;
  double ex;
  double ey;
};

class RadSolver {
 public:
  RadSolver();
  virtual ~RadSolver();
  double getThresholdEnergy() const;
  const TGraphErrors& getBornCrossSection() const;
  const TGraphErrors& getMeasuredCrossSection() const;
  const TMatrixT<double>& getIntegralOeratorMatrix() const;
  const TMatrixT<double>& getInverseErrorMatrix() const;
  bool isThresholdSEnabled() const;
  bool isStartSEnabled() const;
  void solve();
  void setThresholdEnergy(double);
  void setStartPointEnergy(double);
  void setStartPoint(double);
  void setMeasuredCrossSection(const TGraphErrors*);
  void setLeftSideOfBornCrossSection(const TF1*);
  void disableThreshold();
  void enableThreshold();
  void disableStartPoint();
  void enableStartPoint();
  void save(const std::string&);

 private:
  Eigen::MatrixXd getEqMatrix() const;
  double getX(int, int) const;
  static std::pair<double, double> coeffs(double, double, double);
  void check();
  bool _threshold;
  bool _start_point;
  double _threshold_energy;
  double _start_point_enrgy;
  TGraphErrors _measured_cs;
  TGraphErrors _born_cs;
  TF1* _left_side_bcs;
  TMatrixT<double> _inverse_error_matrix;
  TMatrixT<double> _integral_operator_matrix;
  std::vector<RightPart> _measured_cs_data;
};
#endif
