# BICCN
The goal of the BICCN Github repository is to 


The datasets have been made publicly available at LINKforNEMO.
The algorithm used for data integration is online iNMF, detailed in the publication ####. 
To download online iNMF, 
A README file that includes:
1. System requirements
All software dependencies and operating systems (including version numbers)
Versions the software has been tested on
Any required non-standard hardware
2. Installation guide
Instructions
Typical install time on a "normal" desktop computer
3. Demo

Instructions to run on data
Expected output
Expected run time for demo on a "normal" desktop computer
4. Instructions for use
How to run the software on your data
(OPTIONAL) Reproduction instructions




Repository to elicite reproducable workflows for the BICCN project



Compiled standalone software and/or source code
A small (simulated or real) dataset to demo the software/code
A README file that includes:
1. System requirements
All software dependencies and operating systems (including version numbers)
Versions the software has been tested on
Any required non-standard hardware
2. Installation guide
Instructions
Typical install time on a "normal" desktop computer
3. Demo
Instructions to run on data
Expected output
Expected run time for demo on a "normal" desktop computer
4. Instructions for use
How to run the software on your data
(OPTIONAL) Reproduction instructions



## System Requirements

### Hardware requirements
The `rliger` package requires only a standard computer with enough RAM to support the in-memory operations. For minimal performance, please make sure that the computer has at least about 2 GB of RAM. For optimal performance, we recommend a computer with the following specs:

* RAM: 16+ GB
* CPU: 4+ cores, 2.3 GHz/core

### Software requirements

The package development version is tested on *Linux* operating systems and *Mac OSX*.

* Linux: CentOS 7, Manjaro 5.3.18
* Mac OSX: Mojave (10.14.1), Catalina (10.15.2)

The `rliger` package should be compatible with Windows, Mac, and Linux operating systems.

Before setting up the `rliger` package, users should have R version 3.4.0 or higher, and several packages set up from CRAN and other repositories. The user can check the dependencies in `DESCRIPTION`.

## Installation

LIGER is written in R and is also available on the Comprehensive R Archive Network (CRAN). Note that the package name is `rliger` to avoid a naming conflict with an unrelated package. To install the version on CRAN, follow these instructions:

1. Install [R](https://www.r-project.org/)  (>= 3.6)
2. Install [Rstudio](https://posit.co/download/rstudio-desktop/) (recommended)
3. Type the following R command:
```
install.packages('rliger')
```
To install the latest development version directly from GitHub, type the following commands instead of step 3:
```
install.packages('devtools')
library(devtools)
install_github('welch-lab/liger')
```
Note that the GitHub version requires installing from source, which may involve additional installation steps on MacOS (see below).

### Additional Steps for Installing LIGER from Source (recommended before step 3)
Installation from CRAN is easy because pre-compiled binaries are available for Windows and MacOS. However, a few additional steps are required to install from source on MacOS/Windows (e.g. Install RcppArmadillo).
(MacOS) Installing RcppArmadillo on R>=3.4 requires Clang >= 4 and gfortran-6.1. For newer versions of R (R>=3.5), it's recommended to follow the instructions in this [post](https://thecoatlessprofessor.com/programming/r-compiler-tools-for-rcpp-on-macos/). Follow the instructions below if you have R version 3.4.0-3.4.4.

1. Install gfortran as suggested [here](https://gcc.gnu.org/wiki/GFortranBinaries)
2. Download clang4 from this [page](https://mac.R-project.org/libs/clang-4.0.0-darwin15.6-Release.tar.gz)
3. Uncompress the resulting zip file and type into Terminal (`sudo` if needed): 
```
mv /path/to/clang4/ /usr/local/ 
```
4. Create `.R/Makevars` file containing following:
```
# The following statements are required to use the clang4 binary
CC=/usr/local/clang4/bin/clang
CXX=/usr/local/clang4/bin/clang++
CXX11=/usr/local/clang4/bin/clang++
CXX14=/usr/local/clang4/bin/clang++
CXX17=/usr/local/clang4/bin/clang++
CXX1X=/usr/local/clang4/bin/clang++
LDFLAGS=-L/usr/local/clang4/lib
```
For example, use the following Terminal commands:
```
cd ~
mkdir .R
cd .R 
nano Makevars
``` 
Paste in the required text above and save with `Ctrl-X`.

### Additional Installation Steps for Online Learning using LIGER
The HDF5 library is required for implementing online learning in LIGER on data files in HDF5 format. It can be installed via one of the following commands:

| System                                    | Command
|:------------------------------------------|:---------------------------------|
|**OS X (using Homebrew or Conda)**                  | `brew install hdf5` or `conda install -c anaconda hdf5`
|**Debian-based systems (including Ubuntu)**| `sudo apt-get install libhdf5-dev` 
|**Systems supporting yum and RPMs**        | `sudo yum install hdf5-devel`

For Windows, the latest HDF5 1.12.0 is available at https://www.hdfgroup.org/downloads/hdf5/.
