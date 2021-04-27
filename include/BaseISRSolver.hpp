#ifndef _BASE_ISRSOLVER_HPP_
#define _BASE_ISRSOLVER_HPP_

#include <Eigen/Dense>
#include <string>
#include <memory>
#include <TGraphErrors.h>
#include <TEfficiency.h>
#include "ISRSolverStructs.hpp"

/**
 * The exception that is thrown when the detection efficiency function has the
 * wrong number of arguments
 */
typedef struct : std::exception {
  const char* what() const noexcept {
    return "[!] Wrong efficiency dimension.\n";
  }
} EfficiencyDimensionException;

/**
 * Solver base class
 */
class BaseISRSolver {
 public:
  /**
   * Constructor
   * @param inputPath an input path to .root file that contains visible
   * cross section in a form of TGraphErrors object and detection
   * efficiency (if needed) in a form of 1D or 2D TEfficiency object
   * @param inputOpts an input options that contain a name of the
   *  detection efficiency TEfficiency object, a name of the visible
   * cross section TGraphErrors object, a threshold energy
   * @see InputOptions
   */
  BaseISRSolver(const std::string& inputPath, const InputOptions& inputOpts);
  /**
   * Constructor
   * @param vcsGraph a visible cross section in a form of TGraphErrors
   * @param thresholdEnergy a threshold energy (GeV)
   */
  BaseISRSolver(TGraphErrors* vcsGraph, double thresholdEnergy);
  /**
   * Constructor
   * @param vcsGraph a visible cross section in a form of TGraphErrors
   * @param eff a detection efficiency in a form of TEfficiency
   * @param thresholdEnergy a threshold energy (GeV)
   */
  BaseISRSolver(TGraphErrors* vcsGraph,
                TEfficiency* eff,
                double thresholdEnergy);
  /**
   * Copy constructor
   */
  BaseISRSolver(const BaseISRSolver&);
  /**
   * Threshold energy getter
   */
  double getThresholdEnergy() const;
  /**
   * Minimum energy getter
   */
  double getMinEnergy() const;
  /**
   * Maximum energy getter
   */
  double getMaxEnergy() const;
  /**
   * Destructor
   */
  virtual ~BaseISRSolver();
  /**
   * Running the algorithm for finding a solution
   */
  virtual void solve() = 0;
  /**
   * Saving results
   * @param outputPath a path to the .root file where the results
   * are saved
   * @param outputOpts an output options that contain a name
   * of the visible cross section TGraphErrors object
   * (in output file) and a name of the visible cross section
   * TGraphErrors name
   */
  virtual void save(const std::string& outputPath,
                    const OutputOptions& outputOpts) = 0;
  /**
   * Numerical solution (Born cross section) getter
   */
  const Eigen::VectorXd& bcs() const;
  /**
   * Center-of-mass energy const getter
   */
  const Eigen::VectorXd& ecm() const;
  /**
   * Center-of-mass energy error getter
   */
  const Eigen::VectorXd& ecmErr() const;
  /**
   * Visible cross section const getter
   */
  const Eigen::VectorXd& vcs() const;
  /**
   * Visible cross section error const getter
   */
  const Eigen::VectorXd& vcsErr() const;
  /**
   * Detection efficiency function getter
   */
  const std::function<double(double, double)>& efficiency() const;
  /**
   * This method is used to check energy spread mode. If this method
   * returns true than energy spread is enabled
   * (center-of-mass energy error is used by solver in this case),
   * energy spread mode is not enabled otherwise.
   */
  bool isEnergySpreadEnabled() const;
  /**
   * Enable energy spread mode.
   */
  void enableEnergySpread();
  /**
   * Disable energy spread mode
   */
  void disableEnergySpread();
  /**
   * This method is used to reset visible cross section.
   * @param vecVCS a vector of visible cross section values
   * at each center-of-mass energy point
   */
  void resetVisibleCS(const Eigen::VectorXd& vecVCS);
  /**
   * This method is used to reset visible cross section errors.
   * @param vecVCSErr a vector of visible cross section errors
   at each center-of-mass energy point.
   */
  void resetVisibleCSErrors(const Eigen::VectorXd& vecVCSErr);
  /**
   * This method is used to reset center-of-mass energy errors.
   * @param vecECMErr a vector of center-of-mass energy errors.
   */
  void resetECMErrors(const Eigen::VectorXd& vecECMErr);

 protected:
  /**
   * This method is used to get number of center-of-mass
   * energy points
   */
  std::size_t _getN() const;
  /**
   * Threshold energy const getter
   */
  double _energyThreshold() const;
  /**
   * Const getter of a threshold energy square
   */
  double _sThreshold() const;
  /**
   * Const geter of a center-of-mass energy square.
   * @param index an index of a center-of-mass energy point.
   */
  double _s(std::size_t index) const;
  /**
   * Center-of-mass energy const const getter.
   * @param index an index of a center-of-mass energy point.
   */
  double _ecm(std::size_t index) const;
  /**
   * Center-of-mass energy error const getter.
   * @param index an index of the center-of-mass energy point.
   */
  double _ecmErr(std::size_t index) const;
  /**
   * Center-of-mass energy non const getter
   */
  Eigen::VectorXd& _ecm();
  /**
   * Center-of-mass energy error non const getter
   */
  Eigen::VectorXd& _ecmErr();
  /**
   * Visible cross section const getter
   */
  const Eigen::VectorXd& _vcs() const;
  /**
   * Visible cross section non const getter
   */
  Eigen::VectorXd& _vcs();
  /**
   * Visible cross section error const getter
   */
  const Eigen::VectorXd& _vcsErr() const;
  /**
   * Visible cross section error non const getter
   */
  Eigen::VectorXd& _vcsErr();
  /**
   * Born cross section non const getter
   */
  Eigen::VectorXd& _bcs();
  /**
   * Visible cross section error getter in a form of covariance matrix
   */
  Eigen::MatrixXd _vcsInvCovMatrix() const;
  /**
   * Initialize visible cross section data
   * @param vcsGraph a visible cross section in a form of TGraphErrors
   */
  void setupVCS(TGraphErrors* vcsGraph);
  /**
   * Initialize a detection efficiency
   */
  void _setupEfficiency() noexcept(false);

 private:
  /**
   * center-of-mass energy spread
   */
  bool _energySpread;
  /**
   * threshold energy
   */
  double _energyT;
  /**
   * number of center-of-mass energy points
   */
  std::size_t _n;
  /**
   * visible cross section data
   */
  CSVecData _visibleCSData;
  /**
   * detection efficiency function
   */
  std::function<double(double, double)> _efficiency;
  /**
   * detection efficiency
   */
  std::shared_ptr<TEfficiency> _tefficiency;
  /**
   * numerical solution (Born cross section)
   */
  Eigen::VectorXd _bornCS;
};

#endif
