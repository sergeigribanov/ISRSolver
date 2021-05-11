#ifndef _ISRSOLVER_TSVD_HPP_
#define _ISRSOLVER_TSVD_HPP_

#include "ISRSolverSLE.hpp"

/**
 * Solver that solving integral equation using truncated
 singular value decomposition (TSVD)
 */
class ISRSolverTSVD : public ISRSolverSLE {
 public:

  ISRSolverTSVD(std::size_t numberOfPoints,
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
   * @param upperIndex an upper index of TSVD
  */
  ISRSolverTSVD(const std::string& inputPath,
                const InputOptions& inputOpts,
                int upperIndex);
  /**
   * Constructor
   * @param vcsGraph a visible cross section in a form of TGraphErrors
   * @param thresholdEnergy a threshold energy
   * @param upperIndex an upper index of TSVD
   */
  ISRSolverTSVD(TGraphErrors* vcsGraph,
                double thresholdEnergy,
                int upperIndex);
  /**
   * Constructor
   * @parm vcsGraph a visible cross section in a from of TGraphErrors
   * @parm eff a detection efficiency in a form of TEfficiency
   * @param upperIndex an upper index of TSVD
   */
  ISRSolverTSVD(TGraphErrors* vcsGraph,
                TEfficiency* eff,
                double thresholdEnergy,
                int upperIndex);
  /**
   * Copy constructor
   */
  ISRSolverTSVD(const ISRSolverTSVD&);
  /**
   * Destructor
   */
  virtual ~ISRSolverTSVD();
  /**
   * Running the algorithm for finding a solution
   */
  virtual void solve() override;
  /**
   * Setter for the upper TSVD index
   * @param upperIndex an upper TSVD index
   */
  void setUpperTSVDIndex(int upperIndex);
  /**
   * Getter for the upper TSVD index
   */
  int getUpperTSVDIndex() const;
  /**
   * This method enables mode in which just upper TSVD index
   is enabled
   */
  void enableKeepOne();
  /**
   * This method disables mode in which just upper TSVD index
   is enabled
   */
  void disableKeepOne();
 private:
  /**
   * Upper TSVD index
   */
  int _upperTSVDIndex;
  //int _truncIndexUpperLimit;
  /**
   * Keep one mode
   * Just upper TSVD index is enabled if _keepOne=true
   */
  bool _keepOne;
  /**
   * !!! TO DO
   */
  Eigen::MatrixXd _mU;
  Eigen::MatrixXd _mV;
  Eigen::VectorXd _mSing;
};

#endif
