# Code for paper "Regularized Principal Spline Functions to Mitigate Spatial Confounding"
This repository contains scripts, files, and figures related to the simulation study and the real data application described in the paper "Regularized Principal Spline Functions to Mitigate Spatial Confounding".

## Contents of the repository
* **simulations.r** is the main script for running all the simulations using the different approaches compared in the paper.

* **folder "func"** contains the main functions for the simulations:
  - **data_generate.r** contains the function to generate the data for the simulations study.
  - **setting.xslx** contains data table with parameters used to generate data. Each row correspond to a configuration.
  - **sim_functions.r** contains the functions to do the simulations using our approach.
  - **sim_spectralAdj.r** contains functions to run simulations using the Spectral Adjustment approach by Guan et al. (2023).
  - **misc.r** contains auxiliary functions.
 
* **folder "results"** contains some R files that are produced by the simulations. These have been used to create the figures.

* **folder "plots"** contains all the figures in the paper.

* **plot results.r** produces figures related to the simulation study.

* **plot results_suppMat.r** produces figures in the Supplementary Web Material. 

* **ozone.r** reproduces the real data application.

* **SSconf_1.0.tar.gz** is a package containing the function `spsl()` to implement our approach. This is written in C++ using the Armadillo library, and is embedded in an R package to allow for parallelization.

#### More details on these files can be found the file "code-details.html"

## First install all required R packages
  source("install_packages.r")
