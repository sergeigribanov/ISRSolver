#ifndef _ITER_ISR_INTERP_SOLVER_HPP_
#define _ITER_ISR_INTERP_SOLVER_HPP_
#include "ISRSolverSLE.hpp"
#include "Interpolator.hpp"

/**
 * Solver that solving integral equation using iterative method
 */
class IterISRInterpSolver : public ISRSolverSLE {
 public:

  IterISRInterpSolver(std::size_t numberOfPoints,
                      double* energy, double* visibleCS,
                      double* energyErr, double* visibleCSErr,
                      double thresholdEnergy,
                      const std::function<double(double, double)>&
                      efficiency);
  /**
   * Constructor
   * @param inputPath an input path to .root file that contains visible
   cross section in a form of TGraphErrors object and detection
   efficiency (if needed) in a form of 1D or 2D TEfficiency object
   * @param inputOpts an input options that contain a name of the
   detection efficiency TEfficiency object, a name of the visible
   cross section TGraphErrors object, a threshold energy
   * @see InputOptions
   */
  IterISRInterpSolver(const std::string& inputPath,
                      const InputOptions& inputOpts);
  /**
   * Constructor
   * @param vcsGraph a visible cross section in a form of TGraphErrors
   * @param thresholdEnergy a threshold energy
   */
  IterISRInterpSolver(TGraphErrors* vcsGraph,
                      double thresholdEnergy);
  /**
   * Constructor
   * @parm vcsGraph a visible cross section in a from of TGraphErrors
   * @parm eff a detection efficiency in a form of TEfficiency
   */
  IterISRInterpSolver(TGraphErrors* vcsGraph,
                      TEfficiency* eff,
                      double thresholdEnergy);
  /**
   * Copy constructor
   */
  IterISRInterpSolver(const IterISRInterpSolver&);
  /**
   * Destructor
   */
  virtual ~IterISRInterpSolver();
  /**
   * Running the algorithm for finding a solution
   */
  virtual void solve() override;
  /**
   * Setter for a number of iterations
   * @param nIter a number of iterations
   */
  void setNumOfIters(std::size_t nIter);
  /**
   * Getter for a number of iterations
   */
  std::size_t getNumOfIters() const;
  /**
   * Saving results
   * @param outputPath a path to the .root file where the results
   are saved
   * @param outputOpts an output options that contain a name
   of the visible cross section TGraphErrors object
   (in output file) and a name of the visible cross section
   TGraphErrors name
  */
  virtual void save(const std::string&, const OutputOptions&) override;

 private:
  /**
   * Number of iterations
   */
  std::size_t _nIter;
  /**
   * Radiative correction (delta)
   */
  Eigen::VectorXd _radcorr;
};

#endif
