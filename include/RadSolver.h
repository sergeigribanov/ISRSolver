#ifndef __RAD_SOLVER_H__
#define __RAD_SOLVER_H__
#include <vector>
#include "TGraphErrors.h"
#include "TMatrixT.h"

struct RightPart {
  double s_;
  double y_;
  double ex_;
  double ey_;
};


class RadSolver {
 public:
  RadSolver();
  virtual ~RadSolver();
  double getThresholdS() const;
  const TGraphErrors& getBornCrossSection() const;
  const TGraphErrors& getMeasuredCrossSection() const;
  const TMatrixT<double>& getIntegralOeratorMatrix() const;
  const TMatrixT<double>& getInverseErrorMatrix() const;
  void solve();
  void setThresholdS(double);
  void setMeasuredCrossSection(TGraphErrors*);
  void save(const std::string&);
private:
  double _s_threshold;
  TGraphErrors _measured_cs;
  TGraphErrors _born_cs;
  TMatrixT<double> _inverse_error_matrix;
  TMatrixT<double> _integral_operator_matrix;
  TMatrixT<double> getEqMatrix () const;
  double getX(int, int) const;
  std::vector<RightPart> _measured_cs_data;
};


#endif
