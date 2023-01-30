#include <algorithm>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include "BaseISRSolver2.hpp"

double* extractECMPointer(BaseISRSolver2* solver) {
  return solver->_visibleCSData.cmEnergy.data();
}

double* extractECMErrPointer(BaseISRSolver2* solver) {
  return solver->_visibleCSData.cmEnergyError.data();
}

double* extractVCSPointer(BaseISRSolver2* solver) {
  return solver->_visibleCSData.cs.data();
}

double* extractVCSErrPointer(BaseISRSolver2* solver) {
  return solver->_visibleCSData.csError.data();
}

double* extractBCSPointer(BaseISRSolver2* solver) {
  return solver->_bornCS.data();
}

/**
 * Initialize visible cross section data
 * @param vcsGraph a visible cross section in a form of TGraphErrors
 */
void BaseISRSolver2::setupVCS(TGraphErrors* vcsGraph) {
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

BaseISRSolver2::BaseISRSolver2(std::size_t numberOfPoints,
                               double* energy, double* visibleCS,
                               double* energyErr, double* visibleCSErr,
                               double thresholdEnergy,
                               const std::function<double(double, std::size_t)>&
                               efficiency) :
    _energyT(thresholdEnergy),
    _n(numberOfPoints),
    _efficiency(efficiency) {
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
BaseISRSolver2::BaseISRSolver2(TGraphErrors* vcsGraph,
                             double thresholdEnergy)
    : _energySpread(false),
      _energyT(thresholdEnergy),
      _efficiency([](double, std::size_t) {return 1.;}) {
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
BaseISRSolver2::BaseISRSolver2(TGraphErrors* vcsGraph,
                               const std::vector<TH1D*>& eff,
                               double thresholdEnergy) :
    BaseISRSolver2(vcsGraph, thresholdEnergy) {
  _effHists = eff;
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
BaseISRSolver2::BaseISRSolver2(const std::string& inputPath,
                             const InputOptions& inputOpts) :
    _energySpread(false),
    _energyT(inputOpts.thresholdEnergy),
    _efficiency([](double, std::size_t) {return 1.;}) {
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
  _effName = inputOpts.efficiencyName;
  fl->Close();
  delete fl;
}

void BaseISRSolver2::push_back_eff(const std::string& path) {
  auto fl = TFile::Open(path.c_str(), "read");
  auto hist = dynamic_cast<TH1D*>(fl->Get(_effName.c_str())->Clone());
  hist->SetDirectory(0);
  _effHists.push_back(hist);
  fl->Close();
  delete fl;
}

/**
 * Copy constructor
 */
BaseISRSolver2::BaseISRSolver2(const BaseISRSolver2& solver) :
  _energySpread(solver._energySpread),
  _energyT(solver._energyT), _n(solver._n),
  _visibleCSData(solver._visibleCSData),
  _efficiency(solver._efficiency),
  _effHists(solver._effHists),
  _bornCS(solver._bornCS) {}

/**
 * Destructor
 */
BaseISRSolver2::~BaseISRSolver2() {
  for (auto el : _effHists) {
    if (el) {
      delete el;
    }
  }
}

/**
 * Initialize a detection efficiency
 */
void BaseISRSolver2::_setupEfficiency() {
  _efficiency = [this](double x, std::size_t index) {
    int bin = (this->_effHists)[index]->FindBin(x);
    return (this->_effHists)[index]->GetBinContent(bin);
  };
}

std::size_t BaseISRSolver2::getN() const { return _n; }

/**
 * This method is used to get number of center-of-mass
 energy points
*/
std::size_t BaseISRSolver2::_getN() const { return _n; }

/**
 * Threshold energy const getter
 */
double BaseISRSolver2::_energyThreshold() const { return _energyT; }

/**
 * Const getter of a threshold energy square
 */
double BaseISRSolver2::_sThreshold() const { return _energyT * _energyT; }

/**
 * Const geter of a center-of-mass energy square.
 * @param index an index of a center-of-mass energy point.
 */
double BaseISRSolver2::_s(std::size_t i) const {
  return _visibleCSData.cmEnergy(i) * _visibleCSData.cmEnergy(i);
}

/**
 * Center-of-mass energy const const getter.
 * @param index an index of a center-of-mass energy point.
 */
double BaseISRSolver2::_ecm(std::size_t i) const {
  return _visibleCSData.cmEnergy(i);
}

/**
 * Center-of-mass energy const getter
 */
const Eigen::VectorXd& BaseISRSolver2::ecm() const {
  return _visibleCSData.cmEnergy;
}

/**
 * Center-of-mass energy non const getter
 */
Eigen::VectorXd& BaseISRSolver2::_ecm() { return _visibleCSData.cmEnergy; }

/**
 * Center-of-mass energy error getter
 */
const Eigen::VectorXd& BaseISRSolver2::ecmErr() const {
  return _visibleCSData.cmEnergyError;
}

/**
 * Center-of-mass energy error non const getter
 */
Eigen::VectorXd& BaseISRSolver2::_ecmErr() {
  return _visibleCSData.cmEnergyError;
}

/**
 * Center-of-mass energy error const getter.
 * @param index an index of the center-of-mass energy point.
 */
double BaseISRSolver2::_ecmErr(std::size_t i) const {
  return  _visibleCSData.cmEnergyError(i);
}

/**
 * Visible cross section const getter
 */
const Eigen::VectorXd& BaseISRSolver2::_vcs() const {
  return _visibleCSData.cs;
}

/**
 * Visible cross section non const getter
 */
Eigen::VectorXd& BaseISRSolver2::_vcs() { return _visibleCSData.cs; }

/**
 * Visible cross section error const getter
 */
const Eigen::VectorXd& BaseISRSolver2::_vcsErr() const {
  return _visibleCSData.csError;
}

/**
 * Visible cross section error non const getter
 */
Eigen::VectorXd& BaseISRSolver2::_vcsErr() { return _visibleCSData.csError; }

/**
 * Visible cross section error getter in a form of covariance matrix
 */
Eigen::MatrixXd BaseISRSolver2::_vcsInvCovMatrix() const {
  return _visibleCSData.csError.array().pow(-2.).matrix().asDiagonal();
}

/**
 * Numerical solution (Born cross section) getter
 */
const Eigen::VectorXd& BaseISRSolver2::bcs() const { return _bornCS; }

/**
 * Born cross section non const getter
 */
Eigen::VectorXd& BaseISRSolver2::_bcs() { return _bornCS; }

/**
 * Detection efficiency function getter
 */
const std::function<double(double, std::size_t)>& BaseISRSolver2::efficiency() const {
  return _efficiency;
}

/**
 * Threshold energy getter
 */
double BaseISRSolver2::getThresholdEnergy() const {
  return _energyT;
}

/**
 * Minimum energy getter
 */
double BaseISRSolver2::getMinEnergy() const {
  return _visibleCSData.cmEnergy(0);
}

/**
 * Maximum energy getter
 */
double BaseISRSolver2::getMaxEnergy() const {
  return _visibleCSData.cmEnergy(_n - 1);
}

/**
 * This method is used to check energy spread mode. If this method
 returns true than energy spread is enabled
 (center-of-mass energy error is used by solver in this case),
 energy spread mode is not enabled otherwise.
*/
bool BaseISRSolver2::isEnergySpreadEnabled() const {
  return _energySpread;
}

/**
 * Enable energy spread mode.
 */
void BaseISRSolver2::enableEnergySpread() {
  _energySpread = true;
}

/**
 * Disable energy spread mode
 */
void BaseISRSolver2::disableEnergySpread() {
  _energySpread = false;
}

/**
 * This method is used to reset visible cross section.
 * @param vecVCS a vector of visible cross section values
 at each center-of-mass energy point
*/
void BaseISRSolver2::resetVisibleCS(const Eigen::VectorXd& vcs) {
  _visibleCSData.cs = vcs;
}

/**
 * This method is used to reset visible cross section errors.
 * @param vecVCSErr a vector of visible cross section errors
 at each center-of-mass energy point.
*/
void BaseISRSolver2::resetVisibleCSErrors(const Eigen::VectorXd& vcsErr) {
  _visibleCSData.csError = vcsErr;
}

/**
 * This method is used to reset center-of-mass energy errors.
 * @param vecECMErr a vector of center-of-mass energy errors.
 */
void BaseISRSolver2::resetECMErrors(const Eigen::VectorXd& ecmErr) {
  _visibleCSData.cmEnergyError = ecmErr;
}

/**
 * Visible cross section const getter
 */
const Eigen::VectorXd& BaseISRSolver2::vcs() const {
  return _visibleCSData.cs;
}

/**
 * Visible cross section error const getter
 */
const Eigen::VectorXd& BaseISRSolver2::vcsErr() const {
  return _visibleCSData.csError;
}

const std::string& BaseISRSolver2::_getEffName() const {
  return _effName;
}
