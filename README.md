## Introduction
### Relationship between Born and visible cross section
Precise measurement of the inclusive cross section of an e+e- annihilation to hadrons is
the main goal of the experiments carried out with the CMD-3 detector at the VEPP-2000 collider. 
This cross section is of interest in connection with the measurement of the fine structure constant 
and the problem of the muon anomalous magnetic moment. The inclusive cross section of an e+e- 
annihilation to hadrons is considered as the sum of the cross sections for exclusive processes 
of e+e- annihilation into different hadronic states. In the experiments carried out with the CMD-3 
detector, the visible cross section <img src="https://render.githubusercontent.com/render/math?math=\large{\sigma_{\rm vis}}">
is measured, while the Born cross section <img src="https://render.githubusercontent.com/render/math?math=\large{\sigma_{\rm Born}}">
is of interest. The visible and Born cross sections are related by the following integral equation:

![equation Kuraev-Fadin](figures/equation1KuraevFadin.png)

where <img src="https://render.githubusercontent.com/render/math?math=\large{\varepsilon(x, s)}"> is a detection efficiency that depends on 
center-of-mass energy <img src="https://render.githubusercontent.com/render/math?math=\large{\sqrt{s}}"> and the energy fraction carried
away by initial state radiation (ISR). Equation (1) were first obtained by E.A. Kuraev amd V.S. Fadin in the work [1].

It should be noted that the center-of-mass energy is known with some accuracy. This accuracy is determined by the spread of particle energies in each beam and the accuracy of energy measurement. Typical standard deviation of the center-of-mass energy in the case of VEPP-2000 is about 1-3 MeV. Thus, instead of equation (1), the relationship between the visible and Born cross sections is given by equation (2):
![equation Kuraev-Fadin](figures/equation2KuraevFadinBlured.png)

where <img src="https://render.githubusercontent.com/render/math?math=\large{\sigma_{E}(s)}"> is the center-of-mass energy standard deviation. Note that for smooth energy dependences of the cross section, the energy spread can be neglected and equation (1) can be used to relate the visible and Born cross sections.

### ISRSolver
The ISRSolver package is a set of utilities for obtaining the Born cross section from the visible cross section data using various methods, as well as a set of utilities for checking the results. These utilities are available in the form of executable files with a set of command line options. Also it is possible to use the ISRSolver package in a custom C++ or python project.

The most common methods for obtaining the Born cross section using the visible cross section data are the iterative calculation of the radiative correction and the visible cross section fit using qquations (1) or (2), where the Born cross section is specified in some model. 

Equations (1) and (2) can be reduced to a system of linear equations. To reach this goal, it is necessary to interpolate the unknown Born cross section and express the interpolation coefficients in terms of the unknown values of the Born cross section at the energy points where the visible cross section is measured. Then the interpolation of the Born cross section must be substituted into equation (1) or (2). After integration, a system of linear equations for the unknown values of the Born cross section is obtained. This system can be solved directly. This method of finding the Born section is referred to below as the **naive method**. However, it is known that the problem of finding a numerical solution to equation 1 or 2 is an ill-posed problem.

### Naive method

![alt text](figures/equation3KuraevFadinSLE.png)

![alt text](figures/equation4KuraevFadinBluredSLE.png)

### Regularization

### Iterative method

### Visible cross section fit

## Usage

### Naive method

### Regularization

### Iterative method


