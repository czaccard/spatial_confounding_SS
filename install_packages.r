packages <- c("tidyverse", "geoR", "mvtnorm", "tictoc", "doParallel", "coda", 
              "Rcpp", "RcppArmadillo", "splines", "akima", "spBayes", "remotes",
              "LaplacesDemon", "GpGp", "latex2exp", "paletteer")

not_installed <- packages[!(packages %in% installed.packages()[ , "Package"])]    # Extract not installed packages
if(length(not_installed)) install.packages(not_installed)

remotes::install_github("const-ae/sparseMatrixStats") # required by mombf

if (!require('mombf', character.only = TRUE)) {
  install.packages("mombf")
}

if (!require('SSconf', character.only = TRUE)) {
  remotes::install_local("SSconf_1.0.tar.gz")
}