#include <iostream>
#include "Integration.hpp"
#include "KuraevFadin.hpp"
#include "CSplineRangeInterpolator2.hpp"

/**
 * Constructor
 * @param rangeIndexMin a minimum index of sub range (a range between
 two neighboring center-of-mass energy points) that used in interpolation
 range
 * @param rangeIndexMax a maximum index of sub range that used in
 interpolation range
 * @param extCMEnergies extended vector of center-of-mass energies (contains
 threshold energy)
*/
CSplineRangeInterpolator2::CSplineRangeInterpolator2(
    int rangeIndexMin, int rangeIndexMax,
    const Eigen::VectorXd& extCMEnergies):
    BaseRangeInterpolator2(rangeIndexMin, rangeIndexMax, extCMEnergies),
    _acc(std::vector<std::shared_ptr<gsl_interp_accel>>(_numberOfSegments)),
    _spline(std::vector<std::shared_ptr<gsl_spline>>(_numberOfSegments)) {
  Eigen::VectorXd yi = Eigen::VectorXd::Zero(extCMEnergies.rows());
  /**
   * Initialize GSL accelerators and spine (one accelerator and
   spline per each center-of-mass energy segment)
   */
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

/**
 * Copy constructor
 */
CSplineRangeInterpolator2::CSplineRangeInterpolator2(
    const CSplineRangeInterpolator2& rinterp):
    BaseRangeInterpolator2(rinterp),
    _acc(rinterp._acc),
    _spline(rinterp._spline) {}

/**
 * Destructor
 */
CSplineRangeInterpolator2::~CSplineRangeInterpolator2() {}

/**
 * Evaluate basis interpolation around center-of-mass energy that
 corresponds to the cross section point with index the csIndex
 * @param csIndex a cross section point index
 * @param energy a center-of-mass energy at which interpolation
 is evaluated
*/
double CSplineRangeInterpolator2::basisEval(
    int csIndex, double energy) const {
  const int index = csIndex - _beginIndex;
  return gsl_spline_eval(_spline[index].get(), energy, _acc[index].get());
}

/**
 * Evaluate basis interpolation derivative around center-of-mass energy that
 corresponds to the cross section point with the csIndex
 * @param csIndex a cross section point index
 * @param energy a center-of-mass energy a which interpolation is evaluated
 */
double CSplineRangeInterpolator2::basisDerivEval(
    int csIndex, double energy) const {
  const int index = csIndex - _beginIndex;
  return gsl_spline_eval_deriv(_spline[index].get(), energy, _acc[index].get());
}

/**
 * Evaluate Kuraev-Fadin convolution with basis interpolation
 * @param energyIndex a center-of-mass energy index
 * @param csIndex a cross section point index
 * @param efficiency a detection efficiency
 */
double CSplineRangeInterpolator2::evalKuraevFadinBasisIntegral(
    int energyIndex, int csIndex,
    const std::function<double(double, std::size_t)>& efficiency) const {
  const int index = csIndex - _beginIndex;
  /**
   * Get center-of-mass energy that corresponds to energyIndex
   */
  const double en = _spline[index].get()->x[energyIndex + 1];
  if (en <= _minEnergy) {
    return 0;
  }
  std::function<double(double)> eff_1d = [efficiency, energyIndex](double x) {
    return efficiency(x, energyIndex);
  };
  std::function<double(double)> fcn =
      [index, this] (double energy) {
        double result =  gsl_spline_eval(this->_spline[index].get(), energy, this->_acc[index].get());
        return result;
      };
  const double x_min = std::max(0., 1 - std::pow(_maxEnergy / en, 2));
  const double x_max = 1 - std::pow(_minEnergy / en, 2);
  /**
   * Evaluate Kuraev-Fadin convolution with basis interpolation
   */
  return convolutionKuraevFadin_1d(en, fcn, x_min, x_max, eff_1d);
}

double CSplineRangeInterpolator2::evalBasisSConvolution(
    int csIndex,
    const std::function<double(double)>& convKernel) const {
  const int index = csIndex - _beginIndex;
  std::function<double(double)> ifcn =
      [index, convKernel, this] (double s) {
        const double en = std::sqrt(s);
        const double result =  gsl_spline_eval(this->_spline[index].get(), en, this->_acc[index].get()) *
                               convKernel(s);
        return result;
      };
  double error;
  return integrate(ifcn, _minEnergy * _minEnergy, _maxEnergy * _maxEnergy, error);
}

/** Evaluate integral of basis interpolation function that correspond to
 * a csIndex-th cross section point
 * @param csIndex a cross section point index
 */
double CSplineRangeInterpolator2::evalIntegralBasis(int csIndex) const {
  const int index = csIndex - _beginIndex;
  std::function<double(double)> fcn =
    [index, this] (double energy) {
      double result =  gsl_spline_eval(this->_spline[index].get(), energy, this->_acc[index].get());
      return result;
      };
  double error;
  /**
   * Integration
   */
  return integrate(fcn, _minEnergy, _maxEnergy, error);
}

double CSplineRangeInterpolator2::evalBasisSConvolution(
    int csIndex,
    const std::function<double(double)>& convKernel,
    double s_min, double s_max) const {
  const int index = csIndex - _beginIndex;
  std::function<double(double)> ifcn =
      [index, convKernel, this] (double s) {
        const double en = std::sqrt(s);
        const double result =  gsl_spline_eval(this->_spline[index].get(), en, this->_acc[index].get()) *
                               convKernel(s);
        return result;
      };
  double error;
  const double s1_min = std::max(_minEnergy * _minEnergy, s_min);
  const double s1_max = std::min(_maxEnergy * _maxEnergy, s_max);
  return integrate(ifcn, s1_min, s1_max, error);
}
