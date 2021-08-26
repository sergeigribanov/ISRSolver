
![](figures/badge_license_gpl3.svg)
![](figures/badge_linux_passed.svg)

`Documentation in progress...`

## Overview

The ISRSolver toolkit is a set of utilities for obtaining the Born cross section using the visible cross section data, as well as a set of utilities for checking the results. The Born cross section can be found using various methods. Some these methods are generally accepted, while other methods, such as the naive method and the Tikhonov regularization method, were first proposed for finding the Born cross section using the visible cross section data in article [https://arxiv.org/abs/2108.07539](https://arxiv.org/abs/2108.07539 "New method for obtaining a Born cross section using visible cross section data from e+e− colliders") and then implemented in the ISRSolver toolkit.

The utilities are available to the user in a form of executable files that can be run with a set of command line options. The ISRSolver can be also used in a custom C++ or Python project.

## Quick start using Docker
1. Install and setup docker and docker-compose
2. Make sure the Docker service is running. You can do this using ```systemctl```, for example:
  ```console
  systemctl status docker.service
  ```
3. Go to the directory where you want to download the ISRSolver source code and run the following console commands:
  ```console
  git clone https://github.com/sergeigribanov/ISRSolver
  cd ISRSolver
  mkdir shared
  docker-compose up -d
  ```
4. After running the previous commands, check that the isrsolver_isrsolver image is in the list of images: 
  ```console
  docker images
  ```
5. Make sure the isrsolver_isrsolver_1 container is running:
  ```console
  docker ps
  ```
6. Find out the ip-address of the container isrsolver_isrsolvr_1:
  ```console
  docker inspect isrsolver_isrsolver_1
  ```
7. Connect to Jupiter Notebook using your internet browser. In order to do this, use the ip-address from the last point and port 8765. For example, if the ip-address is 172.22.0.2, then you should enter the following URL request in the browser: 172.22.0.2:8765.
8. The default password for Jupyter Notebook is  ```qwedcxzas```.
 
**TO-DO: Describe how to run noteboks and how to use the docker container**

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
source /path/to/ISRSolver/installation/bin/env.sh
```
3. Install required python packages and modules:
```console
pip install --user numpy pandas matplotlib seaborn scikit-hep jupyter
```
4. Go to the notebooks directory:
```console
cd /path/to/ISRSolver/source/code/notebooks
```
5. Download archive with data for numerical experiments https://disk.yandex.com/d/XEes97fR2OsAqQ to the ```notebooks``` directory.
6. Decompress the archive:
```console
tar -xf data.tar.gz
```
7. Run Jupyter Notebook:
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
Each inner list in this file describes interpolation type at a certain set of consecutive c.m. energy intervals. ```false``` means that the interpolation is piecewise linear, ```true``` means cubic spline interpolation. The first number in each inner list is the index of the first c.m. energy interval, and the second number is the index of the last interval. Index ```0``` is the index of c.m. energy interval between the threshold energy and the first c.m. energy point
