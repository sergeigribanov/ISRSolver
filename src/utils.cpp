#include <TRandom3.h>
#include "utils.hpp"

Eigen::VectorXd randomDrawVisCS(const Eigen::VectorXd& vcs,
                                const Eigen::VectorXd& vcsErr) {
  Eigen::VectorXd result = vcs;
  for (int i = 0; i < result.size(); ++i) {
    result(i) = gRandom->Gaus(result(i), vcsErr(i));
  }
  return result;
}
