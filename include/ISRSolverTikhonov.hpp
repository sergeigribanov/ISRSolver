#ifndef __ISRSOLVER_TIKHONOV_HPP__
#define __ISRSOLVER_TIKHONOV_HPP__

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
  bool isDerivNorm2RegIsEnabled() const;
  void enableDerivNorm2Regularizator();
  void disableDerivNorm2Regularizator();

 protected:
  /**
   * This method returns derivative operator matrix
   */
  const Eigen::MatrixXd& _getInterpPointWiseDerivativeProjector() const;
  double _evaldKsidLambda(const Eigen::VectorXd&) const;
  double _evald2Ksid2Lambda(const Eigen::VectorXd&,
                           const Eigen::VectorXd&) const;
  void _evalProblemMatrices();
  void _evalInterpPointWiseDerivativeProjector();

 private:
  bool _enabledDerivNorm2Reg;
  double _lambda;
  Eigen::MatrixXd _interpPointWiseDerivativeProjector;
  Eigen::MatrixXd _mF;
  Eigen::MatrixXd _mL;
  Eigen::FullPivLU<Eigen::MatrixXd> _luT;
  Eigen::FullPivLU<Eigen::MatrixXd> _luL;
};

#endif
