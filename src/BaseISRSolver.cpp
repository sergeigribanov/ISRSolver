#include <algorithm>
#include <vector>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include "BaseISRSolver.hpp"

double* extractECMPointer(BaseISRSolver* solver) {
  return solver->_visibleCSData.cmEnergy.data();
}

double* extractECMErrPointer(BaseISRSolver* solver) {
  return solver->_visibleCSData.cmEnergyError.data();
}

double* extractVCSPointer(BaseISRSolver* solver) {
  return solver->_visibleCSData.cs.data();
}

double* extractVCSErrPointer(BaseISRSolver* solver) {
  return solver->_visibleCSData.csError.data();
}

double* extractBCSPointer(BaseISRSolver* solver) {
  return solver->_bornCS.data();
}

/**
 * Initialize visible cross section data
 * @param vcsGraph a visible cross section in a form of TGraphErrors
 */
void BaseISRSolver::setupVCS(TGraphErrors* vcsGraph) {
  /**
   * Initialize number of points
   */
  _n = vcsGraph->GetN();
  std::vector<CSData> visibleCS;
  visibleCS.reserve(_n);
  for (std::size_t i = 0; i < _n; ++i) {
    /**
     * Load visible cross section data from TGraphErrors object
     */
    visibleCS.push_back(
        {.cmEnergy = vcsGraph->GetX()[i],
        .cs = vcsGraph->GetY()[i],
        .cmEnergyError = vcsGraph->GetEX()[i],
        .csError = vcsGraph->GetEY()[i]});
  }
  /**
   * Sorting a visible cross section data in ascending order of
   center-of-mass energy
   */
  std::sort(visibleCS.begin(), visibleCS.end(),
            [](const CSData& x, const CSData& y) { return x.cmEnergy < y.cmEnergy; });
  /**
   * Converting sorted visible cross section data to a vector format
   */
  _visibleCSData = {
    .cmEnergy = Eigen::VectorXd(_n),
    .cmEnergyError = Eigen::VectorXd(_n),
    .cs = Eigen::VectorXd(_n),
    .csError = Eigen::VectorXd(_n)};
  std::transform(visibleCS.begin(), visibleCS.end(),
                 _visibleCSData.cmEnergy.data(),
                 [](const CSData& x) { return x.cmEnergy; });
  std::transform(visibleCS.begin(), visibleCS.end(),
                 _visibleCSData.cmEnergyError.data(),
                 [](const CSData& x) { return x.cmEnergyError; });
  std::transform(visibleCS.begin(), visibleCS.end(),
                 _visibleCSData.cs.data(),
                 [](const CSData& x) { return x.cs; });
  std::transform(visibleCS.begin(), visibleCS.end(),
                 _visibleCSData.csError.data(),
                 [](const CSData& x) { return x.csError; });
}

BaseISRSolver::BaseISRSolver(std::size_t numberOfPoints,
                             double* energy, double* visibleCS,
                             double* energyErr, double* visibleCSErr,
                             double thresholdEnergy,
                             const std::function<double(double, double)>&
                             efficiency) :
    _energyT(thresholdEnergy),
    _n(numberOfPoints),
    _efficiency(efficiency),
    _tefficiency(nullptr) {
  Eigen::VectorXd enV(_n);
  Eigen::VectorXd csV(_n);
  Eigen::VectorXd enErrV(_n);
  Eigen::VectorXd csErrV(_n);
  std::copy(energy, energy + _n,
                 enV.data());
   std::copy(visibleCS, visibleCS + _n,
             csV.data());
   if (energyErr) {
     _energySpread = true;
     std::copy(energyErr, energyErr + _n,
               enErrV.data());
   } else {
     _energySpread = false;
     enErrV = Eigen::VectorXd::Zero(_n);
   }
  std::copy(visibleCSErr, visibleCSErr + _n,
            csErrV.data());
  std::vector<int> index(_n, 0);
  for (std::size_t i = 0; i < _n; ++i) {
    index[i] = i;
  }
  std::sort(index.begin(), index.end(),
            [enV](int i1, int i2)
            { return enV(i1) < enV(i2); });
  _visibleCSData = {
    .cmEnergy = Eigen::VectorXd(_n),
    .cmEnergyError = Eigen::VectorXd(_n),
    .cs = Eigen::VectorXd(_n),
    .csError = Eigen::VectorXd(_n)};
  std::transform(index.begin(), index.end(),
                 _visibleCSData.cmEnergy.data(),
                 [enV](int i) {return enV(i);});
  std::transform(index.begin(), index.end(),
                 _visibleCSData.cmEnergyError.data(),
                 [enErrV](int i) {return enErrV(i);});
  std::transform(index.begin(), index.end(),
                 _visibleCSData.cs.data(),
                 [csV](int i) {return csV(i);});
  std::transform(index.begin(), index.end(),
                 _visibleCSData.csError.data(),
                 [csErrV](int i) {return csErrV(i);});
}

/**
 * Constructor
 * @param vcsGraph a visible cross section in a form of TGraphErrors
 * @param thresholdEnergy a threshold energy (GeV)
 */
BaseISRSolver::BaseISRSolver(TGraphErrors* vcsGraph,
                             double thresholdEnergy)
    : _energySpread(false),
      _energyT(thresholdEnergy),
      _efficiency([](double, double) {return 1.;}),
      _tefficiency(std::shared_ptr<TEfficiency>(nullptr)) {
  /**
   * Initialize a visible cross section data
   */
  setupVCS(vcsGraph);
}

/**
 * Constructor
 * @param vcsGraph a visible cross section in a form of TGraphErrors
 * @param eff a detection efficiency in a form of TEfficiency
 * @param thresholdEnergy a threshold energy (GeV)
 */
BaseISRSolver::BaseISRSolver(TGraphErrors* vcsGraph,
                             TEfficiency* eff,
                             double thresholdEnergy) :
    BaseISRSolver(vcsGraph, thresholdEnergy) {
  _tefficiency = std::shared_ptr<TEfficiency>(
      dynamic_cast<TEfficiency*>(eff->Clone()));
  /**
   * Initialize a detection efficiency
   */
  _setupEfficiency();
}

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
BaseISRSolver::BaseISRSolver(const std::string& inputPath,
                             const InputOptions& inputOpts) :
    _energySpread(false),
    _energyT(inputOpts.thresholdEnergy),
    _efficiency([](double, double) {return 1.;}),
    _tefficiency(std::shared_ptr<TEfficiency>(nullptr)) {
  /**
   * Opening input file that contains a visible cross section and
   detection efficiency
   */
  auto fl = TFile::Open(inputPath.c_str(), "read");
  /**
   * Get pointer to a visible cross section graph
   */
  auto vcsGraph = dynamic_cast<TGraphErrors*>(
      fl->Get(inputOpts.visibleCSGraphName.c_str()));
  /**
   * Initialize a visible cross section
   */
  setupVCS(vcsGraph);
  /**
   * Load a detection efficiency from the input file
   */
  if (inputOpts.efficiencyName.length() > 0) {
    _tefficiency = std::shared_ptr<TEfficiency>(
        dynamic_cast<TEfficiency*>(fl->Get(inputOpts.efficiencyName.c_str())->Clone()));
  }
  fl->Close();
  delete fl;
  /**
   * Initialize a detection efficiency
   */
  if (_tefficiency.get()) {
    _setupEfficiency();
  }
}

/**
 * Copy constructor
 */
BaseISRSolver::BaseISRSolver(const BaseISRSolver& solver) :
  _energySpread(solver._energySpread),
  _energyT(solver._energyT), _n(solver._n),
  _visibleCSData(solver._visibleCSData),
  _efficiency(solver._efficiency),
  _tefficiency(solver._tefficiency),
  _bornCS(solver._bornCS) {}

/**
 * Destructor
 */
BaseISRSolver::~BaseISRSolver() {}

/**
 * Initialize a detection efficiency
 */
void BaseISRSolver::_setupEfficiency() {
  /**
   * Initialize a 1D detection efficiency
   */
  if (_tefficiency.get()->GetDimension() == 1) {
    _efficiency = [this](double, double energy) {
      int bin = this->_tefficiency.get()->FindFixBin(energy);
      return this->_tefficiency.get()->GetEfficiency(bin);
    };
  }
  /**
   * Initialize a 2D detection efficiency
   */
  if (_tefficiency.get()->GetDimension() == 2) {
    _efficiency = [this](double x, double energy) {
      int bin = this->_tefficiency.get()->FindFixBin(x, energy);
      return this->_tefficiency.get()->GetEfficiency(bin);
    };
  }
  /**
   * Throw an exception if a detection efficiency dimension is wrong
   */
  if (_tefficiency.get()->GetDimension() > 2) {
    EfficiencyDimensionException ex;
    throw ex;
  }
}

std::size_t BaseISRSolver::getN() const { return _n; }

/**
 * This method is used to get number of center-of-mass
 energy points
*/
std::size_t BaseISRSolver::_getN() const { return _n; }

/**
 * Threshold energy const getter
 */
double BaseISRSolver::_energyThreshold() const { return _energyT; }

/**
 * Const getter of a threshold energy square
 */
double BaseISRSolver::_sThreshold() const { return _energyT * _energyT; }

/**
 * Const geter of a center-of-mass energy square.
 * @param index an index of a center-of-mass energy point.
 */
double BaseISRSolver::_s(std::size_t i) const {
  return _visibleCSData.cmEnergy(i) * _visibleCSData.cmEnergy(i);
}

/**
 * Center-of-mass energy const const getter.
 * @param index an index of a center-of-mass energy point.
 */
double BaseISRSolver::_ecm(std::size_t i) const {
  return _visibleCSData.cmEnergy(i);
}

/**
 * Center-of-mass energy const getter
 */
const Eigen::VectorXd& BaseISRSolver::ecm() const {
  return _visibleCSData.cmEnergy;
}

/**
 * Center-of-mass energy non const getter
 */
Eigen::VectorXd& BaseISRSolver::_ecm() { return _visibleCSData.cmEnergy; }

/**
 * Center-of-mass energy error getter
 */
const Eigen::VectorXd& BaseISRSolver::ecmErr() const {
  return _visibleCSData.cmEnergyError;
}

/**
 * Center-of-mass energy error non const getter
 */
Eigen::VectorXd& BaseISRSolver::_ecmErr() {
  return _visibleCSData.cmEnergyError;
}

/**
 * Center-of-mass energy error const getter.
 * @param index an index of the center-of-mass energy point.
 */
double BaseISRSolver::_ecmErr(std::size_t i) const {
  return  _visibleCSData.cmEnergyError(i);
}

/**
 * Visible cross section const getter
 */
const Eigen::VectorXd& BaseISRSolver::_vcs() const {
  return _visibleCSData.cs;
}

/**
 * Visible cross section non const getter
 */
Eigen::VectorXd& BaseISRSolver::_vcs() { return _visibleCSData.cs; }

/**
 * Visible cross section error const getter
 */
const Eigen::VectorXd& BaseISRSolver::_vcsErr() const {
  return _visibleCSData.csError;
}

/**
 * Visible cross section error non const getter
 */
Eigen::VectorXd& BaseISRSolver::_vcsErr() { return _visibleCSData.csError; }

/**
 * Visible cross section error getter in a form of covariance matrix
 */
Eigen::MatrixXd BaseISRSolver::_vcsInvCovMatrix() const {
  return _visibleCSData.csError.array().pow(-2.).matrix().asDiagonal();
}

/**
 * Numerical solution (Born cross section) getter
 */
const Eigen::VectorXd& BaseISRSolver::bcs() const { return _bornCS; }

/**
 * Born cross section non const getter
 */
Eigen::VectorXd& BaseISRSolver::_bcs() { return _bornCS; }

/**
 * Detection efficiency function getter
 */
const std::function<double(double, double)>& BaseISRSolver::efficiency() const {
  return _efficiency;
}

/**
 * Threshold energy getter
 */
double BaseISRSolver::getThresholdEnergy() const {
  return _energyT;
}

/**
 * Minimum energy getter
 */
double BaseISRSolver::getMinEnergy() const {
  return _visibleCSData.cmEnergy(0);
}

/**
 * Maximum energy getter
 */
double BaseISRSolver::getMaxEnergy() const {
  return _visibleCSData.cmEnergy(_n - 1);
}

/**
 * This method is used to check energy spread mode. If this method
 returns true than energy spread is enabled
 (center-of-mass energy error is used by solver in this case),
 energy spread mode is not enabled otherwise.
*/
bool BaseISRSolver::isEnergySpreadEnabled() const {
  return _energySpread;
}

/**
 * Enable energy spread mode.
 */
void BaseISRSolver::enableEnergySpread() {
  _energySpread = true;
}

/**
 * Disable energy spread mode
 */
void BaseISRSolver::disableEnergySpread() {
  _energySpread = false;
}

/**
 * This method is used to reset visible cross section.
 * @param vecVCS a vector of visible cross section values
 at each center-of-mass energy point
*/
void BaseISRSolver::resetVisibleCS(const Eigen::VectorXd& vcs) {
  _visibleCSData.cs = vcs;
}

/**
 * This method is used to reset visible cross section errors.
 * @param vecVCSErr a vector of visible cross section errors
 at each center-of-mass energy point.
*/
void BaseISRSolver::resetVisibleCSErrors(const Eigen::VectorXd& vcsErr) {
  _visibleCSData.csError = vcsErr;
}

/**
 * This method is used to reset center-of-mass energy errors.
 * @param vecECMErr a vector of center-of-mass energy errors.
 */
void BaseISRSolver::resetECMErrors(const Eigen::VectorXd& ecmErr) {
  _visibleCSData.cmEnergyError = ecmErr;
}

/**
 * Visible cross section const getter
 */
const Eigen::VectorXd& BaseISRSolver::vcs() const {
  return _visibleCSData.cs;
}

/**
 * Visible cross section error const getter
 */
const Eigen::VectorXd& BaseISRSolver::vcsErr() const {
  return _visibleCSData.csError;
}
