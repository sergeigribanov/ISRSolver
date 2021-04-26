#ifndef _ISRSOLVER_TIKHONOV_HPP_
#define _ISRSOLVER_TIKHONOV_HPP_

#include "ISRSolverSLE.hpp"

/**
 * Solver that solving integral equation using Tikhonov regularization
 */
class ISRSolverTikhonov : public ISRSolverSLE {
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
   * @param lambda a regularization parameter
   */
  ISRSolverTikhonov(const std::string& inputPath,
                    const InputOptions& inputOpts,
                    double lambda = 1.e-2);
  /**
   * Constructor
   * @param vcsGraph a visible cross section in a form of TGraphErrors
   * @param thresholdEnergy a threshold energy (GeV)
   * @param lambda a regularization parameter
   */
  ISRSolverTikhonov(TGraphErrors* vcsGraph,
                    double thresholdEnergy,
                    double lambda = 1.e-2);
  /**
   * Constructor
   * @param vcsGraph a visible cross section in a form of TGraphErrors
   * @param eff a detection efficiency in a form of TEfficiency
   * @param thresholdEnergy a threshold energy (GeV)
   * @param lambda a regularization parameter
   */
  ISRSolverTikhonov(TGraphErrors* vcsGraph,
                    TEfficiency* eff,
                    double thresholdEnergy,
                    double lambda = 1.e-2);
  /**
   * Copy constructor
   */
  ISRSolverTikhonov(const ISRSolverTikhonov&);
  /**
   * Destructor
   */
  virtual ~ISRSolverTikhonov();
  /**
   * Running the algorithm for finding a solution
   */
  virtual void solve() override;
  /**
   * Regularization parameter getter
   */
  double getLambda() const;
  /**
   * Regularization parameter setter
   * @param lambda a regularization parameter
   */
  void setLambda(double lambda);
  /**
   * This method evaluates the visible cross section chi-square
   */
  double evalEqNorm2() const;
  /**
   * This method evaluates L2 norm of numerical solution or its derivative
   */
  double evalSmoothnessConstraintNorm2() const;
  /**
   * This method evaluates L-curve curvature
   */
  double evalLCurveCurvature() const;
  /**
   * This method evaluates L-curve curvature derivative
   */
  double evalLCurveCurvatureDerivative() const;
  /**
   * This method returns true in the case when L2 norm square of
   the numerical solution is used as regularizator. Otherwise this method
   returns false.
   */
  bool isDerivNorm2RegIsEnabled() const;
  /**
   * This method enables regularizator in a form of L2 norm square
   of the numerical solution derivative and disables regularizator in a
   form of L2 norm square of the numerical solution
   */
  void enableDerivNorm2Regularizator();
  /**
   * This method disables regularizator in a form of L2 norm square
   of the numerical solution derivative and enables regularizator in a
   form of L2 norm square of the numerical solution
   */
  void disableDerivNorm2Regularizator();

 protected:
  /**
   * This method returns derivative operator matrix
   */
  const Eigen::MatrixXd& _getInterpPointWiseDerivativeProjector() const;
  /**
   * !!! TO DO
   */
  double _evaldKsidLambda(const Eigen::VectorXd&) const;
  /**
   * !!! TO DO
   */
  double _evald2Ksid2Lambda(const Eigen::VectorXd&,
                           const Eigen::VectorXd&) const;
  /**
   * This method calculates auxiliary matrices arising from regularization
   */
  void _evalProblemMatrices();
  /**
   * This method evaluates derivative operator matrix
   */
  void _evalInterpPointWiseDerivativeProjector();

 private:
  /**
   * This variable is true if regularizator uses numerical solution
   derivative
   */
  bool _enabledDerivNorm2Reg;
  /**
   * Regularization parameter
   */
  double _lambda;
  /**
   * Derivative operator matrix
   */
  Eigen::MatrixXd _interpPointWiseDerivativeProjector;
  /*
   * !!! TO DO
   */
  Eigen::MatrixXd _mF;
  Eigen::MatrixXd _mL;
  Eigen::FullPivLU<Eigen::MatrixXd> _luT;
  Eigen::FullPivLU<Eigen::MatrixXd> _luL;
};

#endif
