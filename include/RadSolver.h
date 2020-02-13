#ifndef __RAD_SOLVER_H__
#define __RAD_SOLVER_H__
#include <vector>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TMatrixT.h>

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
  TMatrixT<double> getEqMatrix () const;
  double getX(int, int) const;
  std::vector<RightPart> _measured_cs_data;
};


#endif
