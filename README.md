
![](figures/badge_license_gpl3.svg)
![](figures/badge_linux_passed.svg)

`Documentation in progress...`

## Overview

The ISRSolver toolkit is a set of utilities for obtaining a Born cross section using visible cross section data, as well as a set of utilities for checking the results. The Born cross section can be found using various methods. Some these methods are generally accepted, while other methods, such as the naive method and the Tikhonov regularization method, were first proposed for finding the Born cross section using the visible cross section data in article [https://arxiv.org/abs/2108.07539](https://arxiv.org/abs/2108.07539 "New method for obtaining a Born cross section using visible cross section data from e+e− colliders") and then implemented in the ISRSolver toolkit.

The utilities are available to the user in the form of executable files that can be run with a set of command line options. The ISRSolver can be also used in a custom C++ or Python project.

## Quick start using Docker
### Installation using a pre-built image 
1. Install and setup docker and docker-compose.
2. Make sure the Docker service is running. You can do this using ```systemctl```, for example:
  ```console
  systemctl status docker.service
  ```
3. Download ```ssgribanov/isrsolver``` image from the official Docker repository:
```console
docker pull ssgribanov/isrsolver
```
4. Check that the ``ssgribanov/isrsolver``` image is actually downloaded:
```console
docker images
```
5. Go to the directory where you want to download the ISRSolver source code and run the following console commands:
  ```console
  git clone https://github.com/sergeigribanov/ISRSolver
  cd ISRSolver
  mkdir shared
  docker-compose up -d
  ```
6. Check that the container ```isrsolver_isrsolver_1``` is actually running.
 
### Manual build
1. Install and setup docker and docker-compose.
2. Make sure the Docker service is running. You can do this using ```systemctl```, for example:
  ```console
  systemctl status docker.service
  ```
3. Go to the directory where you want to download the ISRSolver source code and run the following console commands:
  ```console
  git clone https://github.com/sergeigribanov/ISRSolver
  cd ISRSolver
  mkdir shared
  ```
4. Start image building:
```console
docker build -t ssgribanov/isrsolver:latest .
```
5. After the build is complete, check that the ```ssgribanov/isrsolver``` image was actually built:
  ```console
  docker images
  ```
6. Run container using the following command:
```console
docker-compose up -d
```
7. Check that the container ```isrsolver_isrsolver_1``` is actually running.
  
### Usage (Jupyter Notebooks)
1. Check the container ```IP```:
```console
docker inspect isrsolver_isrsolver_1 | grep IPAddress
````
The container tag may differ from the standard one (```isrsolver_isrsolver_```) if it was launched in a directory other than ```ISRSolver```. In the previous command, you need to use the actual container tag.
2. Connect to Jupiter Notebook using your internet browser. In order to do this, use the ip-address from the last point and port 8765. For example, if the ip-address is ```172.24.0.2```, then you should enter the following URL request in the browser: ```172.24.0.2:8765```.
4. Enter the default password for Jupyter Notebook is  ```qwedcxzas```.
5. After the previous command, two directories will be available in the browser window: ```notebooks``` and ```shared```.
The ```shared``` directory is intended for files that you want to export or import into the container. This directory corresponds to the host directory created at the step 5 of the ***Installation using a pre-built image*** section or at the step 3 of the ***Manual build*** section. The ```notebooks``` directory contains a collection of Jupiter Notebooks that you can open and run.

### Usage (Console)
1. Run the ```/bin/bash``` command inside the existing container:
```console
docker exec -it isrsolver_isrsolver_1 /bin/bash
```
2. Now you able run any ISRSolver executable files, for example:
```console
isrsolver-SLE -t 0.827 -i notebooks/data/gen_visible_cs_etapipi_simple_model_no_energy_spread.root -o shared/test.root
```
Detailed information on running ISRSolver executable files is given below in section ***Usage***

## Installation
1. Make sure that packages [ROOT](https://root.cern "ROOT - Data Analysis Framework") (```C++11```), [GSL](https://www.gnu.org/software/gsl "GSL - GNU Scientific Library"), [Eigen 3](https://eigen.tuxfamily.org/index.php?title=Main_Page "Eigen - C++ template library for linear algebra"), [Boost](https://www.boost.org "Boost - free peer-reviewed portable C++ source libraries"), [NLopt](https://nlopt.readthedocs.io/en/latest "NLopt - free/open-source library for nonlinear optimization"), [nlohmann_json](https://github.com/nlohmann/json "JSON for Modern C++") and [Minuit2 stand-alone](https://github.com/GooFit/Minuit2 "Stand-alone Minuit2") are installed.
2. In the following we assume that you have the following directory tree:
 - ```$HOME/source``` - the source code directory,
 - ```$HOME/build``` - the build directory,
 - ```$HOME/packages``` - the installation directory.

3. Download the ISRSolver source code:
  ```console
  git clone https://github.com/sergeigribanov/ISRSolver $HOME/source/ISRSolver
  ```
4. Create a directory ```$HOME/build/ISRSolver``` and change to it:
  ```console
  mkdir $HOME/build/ISRSolver
  cd $HOME/build/ISRSolver
  ```
5. Setup ```ROOT``` environment.
6. Run the following command:
  ```console
  cmake -DCMAKE_INSTALL_PREFIX=$HOME/packages/ISRSolver $HOME/source/ISRSolver
  ```
7. Please note that cmake sometimes cannot find some packages depending on how they were installed. If an error occurs that cmake cannot find a particular package, you should to run cmake with the appropriate options. **TO-DO: insert these commands and describe options**
8. Build and install:
  ```console
  make -j8
  make install
  ```
## Jupyter Notebooks without Docker
1. Setup ```ROOT``` environment.
2. Setup ```ISRSolver``` environment:
```console
source <ISRSolver installation prefix>/bin/env.sh
```
3. Install required python (Python 3) packages and modules:
```console
pip install --user numpy pandas matplotlib seaborn scikit-hep jupyter
```
4. Go to the notebooks directory:
```console
cd <ISRSolver source code prefix>/notebooks
```
5. Make directory ```../shared```. This directory is used in some examples to save results.
```console
mkdir ../shared
```
6. Download archive with data for numerical experiments https://disk.yandex.com/d/XEes97fR2OsAqQ to the ```notebooks``` directory.
7. Decompress the archive:
```console
tar -xf data.tar.gz
```
8. Run Jupyter Notebook:
```console
jupyter notebook
```

## Usage
Suppose you have a ```.root``` file that stores a visible cross section in the form of ```TGraphErrors``` object and a detection efficiency in the form of 2D ```TEfficiency``` object. Suppose this file is called ```input.root```, the visible cross section graph is called ```vcs``` and the detection efficiency object is called ```efficiency```. In what follows, it is assumed that the vertical error bars on the ```vcs``` graph represent visible cross section uncertainties, and the horizontal error bars represent c.m. energy spread. The detection efficiency is a function of two variables. The first variable is ```x```, while the second variable is a ```c.m. energy```.

### Naive method
1. Setup ```ROOT``` environment:
```console
source <path to ROOT installation>/bin/thisroot.sh
```
2. Setup ```ISRSolver``` environment:
``` console
source <path to ISRSolver installation>/bin/env.sh
```
3. Run the following command:
```console
isrsolver-SLE -t <threshold energy> -g -v vcs -e efficiency -i <path to input.root> -o output.root
```
Option ```-t``` is used to set a threshold energy in GeV, option ```-i``` is used to set a path to the input file and option ```-o``` is used to set an output file path. Option ```-v``` is used to set the name of the visible cross section graph stored in the input file. If this name is ```vcs``` this option can be omitted. Option ```-e``` is used to set the name of the detection efficiency object stored in the input file. If the option ```-e``` is omitted, then it is assumed that the detection is equalt to ```1``` or the visible cross section is already divided by the one-dimensional detection efficiency. Option ```-g``` is used to take into account c.m. energy spread. If this option is omitted, the c.m. energy spread is not taken into account even if the visible cross section graph has horizontal error bars.

By default, the Born section is interpolated using a piecewise linear interpolation. To change the type of interpolation, option ```-r``` can be used to set a path to a ```.json``` file with interpolation settings:
```console
isrsolver-SLE -t <threshold energy> -g -v vcs -e efficiency -i <path to input.root> -o output.root -r interp-settings.json
```
Suppose there are 50 c.m. energy points at which a visible cross-section is measured. Suppose also that we want the interpolation to be piecewise linear on the interval from the threshold energy to the first measurement point, and interpolation using a cubic spline at the intervals between the other measurement points. In this case, the config file ```interp-settings.json``` will be like this
```json
[
[false, 0, 0], 
[true, 1, 49]
]
```
Each inner list in this file describes interpolation type at a certain set of consecutive c.m. energy intervals. ```false``` means that the interpolation is piecewise linear, ```true``` means cubic spline interpolation. The first number in each inner list is the index of the first c.m. energy interval, and the second number is the index of the last interval. Index ```0``` is the index of c.m. energy interval between the threshold energy and the first c.m. energy point.

### Tikhonov regularization
Regularization leads to a biased numerical solution, so the covariance matrix of the Born cross section is incorrect. 
#### Solving using the manual regularization parameter 
1. Setup ```ROOT``` and ```ISRSolver``` environments.
2. Run the following command:
```console
isrsolver-Tikhonov -t <threshold energy> -l 1.0 -i <path to input.root> -o output.root
```
Option ```-l``` is used to set the regularization parameter manually. Options ```-t```, ```-i``` and ```-o``` are the same as in the case of the ```isrsolver-SLE``` utility. Options ```-g```, ```-e```, ```-v``` and ```-r``` are also the same as in the case of the ```isrsolver-SLE``` utility and should be used if needed. By default, the regularization term contains the square of the norm of the numerical solution derivative. To use the square of the norm of the numerical solution as a regularization term, option ```-s``` should be enabled.
#### L-Curve and L-Curve curvature plots
1. Setup ```ROOT``` and ```ISRSolver``` environments.
2. Run the following command:
```console
isrsolver-LCurve-plot -t <threshold energy> -m <min reg. param.> -x <max reg. param.> -n <number of reg. param. points> -i <path to input.root> -o output.root
```
Options ```-m``` and ```-x``` are used to set minimum and maximum values of the regularization parameter, respectively. Option ```-n``` is used to set number of steps in the regularization parameter. Options ```-t```, ```-i``` and ```-o``` are the same as in the case of ```isrsolver-Tikhonov``` utility. Options ```-g```, ```-e```, ```-v```, ```-r``` and ```-s``` are also the same as in the case of ```isrsolver-Tikhonov``` utility and should be used if needed.
#### Solving using the L-Curve criterion
1. Setup ```ROOT``` and ```ISRSolver``` environments.
2. Run the following command:
```console
isrsolver-Tikhonov-LCurve -t <threshold energy> -l <initial reg. param.> -i <path to input.root> -o output.root
```
Option ```-l``` is used to set an initial value of the regularization parameter. It is better to choose this value close to the point of maximum curvature of the L-Curve. Options ```-t```,  ```-i``` and ```-o``` are the same as in the case of the ```isrsolver-Tikhonov``` utility. Options ```-g```, ```-e```, ```-v```, ```-r``` and ```-s``` are also the same as in the case of the ```isrsolver-Tikhonov``` utility and should be used if needed.

### Validation
The numerical solution validation can be performed using the model data. The model data should be prepared as follows:
1. Model dependence of the Born crosss section on c.m. energy is known.
2. The model visible cross section is calculated using the model Born cross section, efficiency, and energy spread at a given c.m. energy points.
3. "Experimental visible cross section" values at points are generated according to the model visible cross section and a given uncertainties of the visible cross section.
4. The "experimental visible cross section" and the detection efficiency are recorded to a file (cross section - ```TGraphErrors```, efficiency - ```TEfficiecny```).
5. The model Born cross section is evaluated at the given c.m. energy points and stored in the form of ```TGraphErrors``` to a file. The model visible cross section is also stored to the same file in the form of ```TGraphErrors```. 
#### Chi-square test
##### Naive method
1. Setup ```ROOT``` and ```ISRSolver``` environments.
2. Run the following command:
```console
isrsolver-SLE-chi2-test -t <threshold energy> -u <path to the file with model Born and visible cross sections> -b <name of the model Born C.S. graph> -c <name of the model visible C.S. graph> -i <path to the file with the "experimental visible C.S." and the detection eff. (if the efficiecny is used)> -v <name of the "experimental visible C.S. graph> -e <name of the detection eff. object (if the efficiency is used)> -n <number of toy Monte Carlo events> -o output.root
```
Option ```-u``` is used to set the path to the file with model Born and visible cross section graphs. Option ```-b``` is used to set the name of the model Born cross section graph in this file, while option ```-c``` is used to set the name of the model visible cross section, which is also stored in this file. Option ```-i``` is used to set the path to the "experimental visible cross section" and the detection efficiency (in case if the detection efficiency is used). Option ```-v``` is used to set the name of the "experimental visible cross section" graph, while option ```-e``` is used to set the name of the detection efficiency object. In the case if the detection efficiency not used, the option ```-e``` should be omitted. Option ```-n``` is used to set a number of toy Monte Carlo evants. The utility creates a chi-square histogram. This histogram is then fitted with the appropriate distribution. To set the initial amplitude of this distribution, option ```-a``` can be used. Options ```-g``` and ```-r``` are the same as in the case of the ```isrsolver-SLE``` utility and should be used if needed.
##### Tikhonov regularization
1. Setup ```ROOT``` and ```ISRSolver``` environments.
2. Run the following command:
```console
isrsolver-Tikhonov-chi2-test -l <reg. param.> -t <threshold energy> -u <path to the file with model Born and visible cross sections> -b <name of the model Born C.S. graph> -c <name of the model visible C.S. graph> -i <path to the file with the "experimental visible C.S." and the detection eff. (if the efficiecny is used)> -v <name of the "experimental visible C.S. graph> -e <name of the detection eff. object (if the efficiency is used)> -n <number of toy Monte Carlo events> -o output.root
```
Option ```-l``` is used to set the regularization parameter. Other options are the same as in the ```isrsolver-SLE-chi2-test``` utility. Options ```-a```, ```-g```, ```-r``` and ```-s``` can be used if needed.
#### Ratio test
##### Naive method
1. Setup ```ROOT``` and ```ISRSolver``` environments.
2. Run the following command:
```console
isrsolver-SLE-ratio-test -t <threshold energy> -n <number of toy Monte Carlo events> -u <path to the file with model Born and visible cross sections> -b <name of the model Born C.S. graph> -c <name of the model visible C.S. graph> -i <path to the file with the "experimental visible C.S." and the detection eff. (if the efficiecny is used)> -v <name of the "experimental visible C.S. graph> -e <name of the detection eff. object (if the efficiency is used)> -o output.root
```
Options ```-g``` and ```-r``` are the same as in the case of the ```isrsolver-SLE``` utility and should be used if needed.
##### Tikhonov regularization
1. Setup ```ROOT``` and ```ISRSolver``` environments.
2. Run the following command:
```console
isrsolver-Tikhonov-ratio-test -l <reg. param.> -t <threshold energy> -n <number of toy Monte Carlo events> -u <path to the file with model Born and visible cross sections> -b <name of the model Born C.S. graph> -c <name of the model visible C.S. graph> -i <path to the file with the "experimental visible C.S." and the detection eff. (if the efficiecny is used)> -v <name of the "experimental visible C.S. graph> -e <name of the detection eff. object (if the efficiency is used)> -o output.root
```
Option ```-l``` is used to set the regularization parameter. Other options are the same as in the ```isrsolver-SLE-chi2-test``` utility. Options ```-g```, ```-r``` and ```-s``` can be used if needed.
