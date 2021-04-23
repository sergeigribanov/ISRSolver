#ifndef __ISRSOLVER_SLE_HPP__
#define __ISRSOLVER_SLE_HPP__

#include <TF1.h>
#include "Interpolator.hpp"
#include "BaseISRSolver.hpp"

/**
 * Solver that solving integral equation using the naive method (without
 regularization). In this case the integral equation is reduced to a system of
 linear equations and then these equations are solved directly
 */
class ISRSolverSLE : public BaseISRSolver {
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
  ISRSolverSLE(const std::string& inputPath,
               const InputOptions& inputOpts);
  /**
   * Constructor
   * @param vcsGraph a visible cross section in a form of TGraphErrors
   * @param thresholdEnergy a threshold energy (GeV)
   */
  ISRSolverSLE(TGraphErrors* vcsGraph,
               double thresholdEnergy);
  /**
   * Constructor
   * @param vcsGraph a visible cross section in a form of TGraphErrors
   * @param eff a detection efficiency in a form of TEfficiency
   * @param thresholdEnergy a threshold energy (GeV)
   */
  ISRSolverSLE(TGraphErrors* vcsGraph,
               TEfficiency* eff,
               double thresholdEnergy);
  /**
   * Copy constructor
   */
  ISRSolverSLE(const ISRSolverSLE&);
  /**
   * Destructor
   */
  virtual ~ISRSolverSLE();
  /**
   * This method returns integral operator matrix
   */
  const Eigen::MatrixXd& getIntegralOperatorMatrix() const;
  /**
   * This method returns covariance matrix of numerical solution
   (Born cross section)
   */
  const Eigen::MatrixXd& getBornCSCovMatrix() const;
  /**
   * Running the algorithm for finding a solution
   */
  virtual void solve() override;
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
                    const OutputOptions& outputOpts) override;
  /**
   * This method applying interpolation settings of numerical solution
   * @param interpRangeSettings an interpolation settings
   */
  void setRangeInterpSettings(const std::vector<std::tuple<bool, int, int>>&
                              interpRangeSettings);
  /**
   * This method applying interpolation settings of numerical solution
   * @param pathToJSON a path to .json file with interpolation settings
   */
  void setRangeInterpSettings(const std::string& pathToJSON);
  /**
   * This method computes integral operator matrix
   */
  void evalEqMatrix();
  /**
   * This method evaluates the value of the interpolation function at
   a certain energy
   * @param y an interpolated points (Y-values)
   * @param energy a center-of-mass energy at which interpolation function is evaluating
   */
  double interpEval(const Eigen::VectorXd& y, double energy) const;
  /**
   * This method evaluates condition number of the integral operator matrix
   */
  double evalConditionNumber() const;
  /**
   * This method prints condition number of (non regularized) integral operator matrix
   */
  void printConditionNumber() const;

 protected:
  /**
   * Integral operator matrix non const getter
   */
  Eigen::MatrixXd& _getIntegralOperatorMatrix();
  /**
   * Covariance matrix non const getter
   */
  Eigen::MatrixXd& _getBornCSCovMatrix();
  /**
   * This method returns center-of-mass energy spread matrix
   */
  Eigen::MatrixXd _energySpreadMatrix() const;
  /**
   * This method returns dot product weights
   */
  const Eigen::RowVectorXd& _getDotProdOp() const;
  /**
   * This method evaluates dot product operator
   */
  void _evalDotProductOperator();
  /**
   * This method creates interpolation function for the numerical
   solution (Born cross section) in a form of TF1
   */
  TF1* _createInterpFunction() const;
  /**
   * Numerical solution (Born cross section) error const getter
   */
  Eigen::VectorXd _bcsErr() const;
  /**
   * Interpolator that interpolates the numerical solution
   */
  Interpolator _interp;
  /**
   * A boolean flag that is true when integral operator matrix is
   prepared and false otherwise
   */
  bool _isEqMatrixPrepared;

 private:
  /**
   * Integral operator matrix
   */
  Eigen::MatrixXd _integralOperatorMatrix;
  /**
   * Covariance matrix of the numerical solution
   */
  Eigen::MatrixXd _covMatrixBornCS;
  /**
   * Dot product weights
   */
  Eigen::RowVectorXd _dotProdOp;
};

#endif
