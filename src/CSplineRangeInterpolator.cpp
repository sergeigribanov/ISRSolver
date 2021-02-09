#include "kuraev_fadin.hpp"
#include "CSplineRangeInterpolator.hpp"

CSplineRangeInterpolator::CSplineRangeInterpolator(
    int rangeIndexMin, int rangeIndexMax,
    const Eigen::VectorXd& extCMEnergies):
    BaseRangeInterpolator(rangeIndexMin, rangeIndexMax, extCMEnergies),
    _acc(std::vector<std::shared_ptr<gsl_interp_accel>>(_numberOfSegments)),
    _spline(std::vector<std::shared_ptr<gsl_spline>>(_numberOfSegments)) {
  Eigen::VectorXd yi = Eigen::VectorXd::Zero(extCMEnergies.rows());
  for (int i = 0; i < _numberOfSegments; ++i) {
    _acc[i] = std::shared_ptr<gsl_interp_accel>(
        gsl_interp_accel_alloc(),
        [](gsl_interp_accel* tacc)
        {gsl_interp_accel_free(tacc);});
    _spline[i] = std::shared_ptr<gsl_spline>(
        gsl_spline_alloc(gsl_interp_cspline, extCMEnergies.rows()),
        [](gsl_spline* tspline)
        {gsl_spline_free(tspline);});
    yi(_beginIndex + i + 1) = 1.;
    gsl_spline_init(_spline[i].get(), extCMEnergies.data(),
                    yi.data(), extCMEnergies.rows());
    yi(_beginIndex + i + 1) = 0.;
  }
}

CSplineRangeInterpolator::CSplineRangeInterpolator(
    const CSplineRangeInterpolator& rinterp):
    BaseRangeInterpolator(rinterp),
    _acc(rinterp._acc),
    _spline(rinterp._spline) {}

CSplineRangeInterpolator::~CSplineRangeInterpolator() {}

double CSplineRangeInterpolator::basisEval(
    int csIndex, double energy) const {
  const int index = csIndex - _beginIndex;
  return gsl_spline_eval(_spline[index].get(), energy, _acc[index].get());
}

double CSplineRangeInterpolator::basisDerivEval(
    int csIndex, double energy) const {
  const int index = csIndex - _beginIndex;
  return gsl_spline_eval_deriv(_spline[index].get(), energy, _acc[index].get());
}

double CSplineRangeInterpolator::evalKuraevFadinBasisIntegral(
    int energyIndex, int csIndex,
    const std::function<double(double, double)>& efficiency) const {
  const int index = csIndex - _beginIndex;
  std::function<double(double)> fcn =
      [index, this] (double energy) {
        double result =  gsl_spline_eval(this->_spline[index].get(), energy, this->_acc[index].get());
        return result;
      };
  const double x_min = 1 - std::pow(_spline[index].get()->x[0] / _minEnergy, 2);
  const double x_max = 1 - std::pow(_spline[index].get()->x[0] / _spline[index].get()->x[energyIndex + 1], 2);
  double integralResult = kuraev_fadin_convolution(_spline[index].get()->x[energyIndex + 1], fcn, x_min, x_max, efficiency);
  return integralResult;
}
