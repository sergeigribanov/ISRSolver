#ifndef __BASE_ISRSOLVER_HPP__
#define __BASE_ISRSOLVER_HPP__

#include <Eigen/Dense>
#include <string>
#include <TEfficiency.h>

#include "ISRSolverStructs.hpp"

/**
 * The exception that is thrown when the detection efficiency function has the
 wrong number of arguments
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
   cross section in a form of TGraphErrors object and detection
   efficiency (if needed) in a form of 1D or 2D TEfficiency object
   * @param inputOpts an input options that contain a name of the
   detection efficiency TEfficiency object, a name of the visible
   cross section TGraphErrors object, a threshold energy
   @see InputOptions
   */
  BaseISRSolver(const std::string& inputPath, const InputOptions& inputOpts);
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
   are saved
   * @param outputOpts an output options that contain a name
   of the visible cross section TGraphErrors object
   (in output file) and a name of the visible cross section
   TGraphErrors name
   */
  virtual void save(const std::string& outputPath,
                    const OutputOptions& outputOpts) = 0;
  const Eigen::VectorXd& bcs() const;
  const Eigen::VectorXd& ecm() const;
  const Eigen::VectorXd& ecmErr() const;
  const Eigen::VectorXd& vcs() const;
  const Eigen::VectorXd& vcsErr() const;
  const std::function<double(double, double)>& efficiency() const;
  bool isEnergySpreadEnabled() const;
  void enableEnergySpread();
  void disableEnergySpread();
  void resetVisibleCS(const Eigen::VectorXd&);
  void resetVisibleCSErrors(const Eigen::VectorXd&);
  void resetECMErrors(const Eigen::VectorXd&);

 protected:
  std::size_t _getN() const;
  double _energyThreshold() const;
  double _sThreshold() const;
  double _s(std::size_t) const;
  double _ecm(std::size_t) const;
  double _ecmErr(std::size_t) const;
  Eigen::VectorXd& _ecm();
  Eigen::VectorXd& _ecmErr();
  const Eigen::VectorXd& _vcs() const;
  Eigen::VectorXd& _vcs();
  const Eigen::VectorXd& _vcsErr() const;
  Eigen::VectorXd& _vcsErr();
  Eigen::VectorXd& _bcs();
  Eigen::MatrixXd _vcsInvCovMatrix() const;
  void _setupEfficiency() noexcept(false);

 private:
  bool _energySpread;
  double _energyT;
  std::size_t _n;
  CSVecData _visibleCSData;
  std::function<double(double, double)> _efficiency;
  TEfficiency* _tefficiency;
  Eigen::VectorXd _bornCS;
};

#endif
