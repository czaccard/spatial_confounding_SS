packages <- c("tidyverse", "geoR", "mvtnorm", "tictoc", "parallel", "doParallel", "coda", 
              "Rcpp", "RcppArmadillo", "splines", "akima", "mombf", "spBayes", "remotes",
              "LaplacesDemon", "GpGp", "latex2exp")


not_installed <- packages[!(packages %in% installed.packages()[ , "Package"])]    # Extract not installed packages
if(length(not_installed)) install.packages(not_installed)

if (!require('SSconf', character.only = TRUE)) {
  remotes::install_local("SSconf_1.0.tar.gz")
}