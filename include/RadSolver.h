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
  RadSolver(TGraphErrors*, double);
  ~RadSolver();
  TGraphErrors* getBornCS(TMatrixT<double>&,
			  TMatrixT<double>&) const;
  TGraphErrors* visible_cs;
 private:
  TMatrixT<double> getEqMatrix () const;
  double getX(int, int) const;
  std::vector<RightPart> rightPart_;
};


#endif
