rm(list=ls())
outpath = "results/"
inpath = "func/"

library(tidyverse)
library(geoR)
library(mvtnorm)
library(tictoc)
library(parallel)
library(doParallel)
library(coda)
library(Rcpp)
library(splines)
library(akima)
library(SSconf)
library(mombf)
library(spBayes)
library(mgcv)
library(metR)
library(paletteer)

source(paste0(inpath, "data_generate.r"), local = T)
source(paste0(inpath, "sim_functions.r"), local = T)
source(paste0(inpath, "misc.r"), local = T)
source(paste0(inpath, "sim_spectralAdj.r"), local = T)

beta.real = c(1,2) %>% as.matrix(ncol=1) # vector giving intercept coefficient and beta_x
ind = 2 # 2: indicates the exposure coefficient in a design matrix like cbind(1,X)

M = 2^6 # Length of one side of grid 
griddf = expand.grid(lngcoords = seq(0, 1, l=M+1)[-(M+1)],
                     latcoords = seq(0, 1, l=M+1)[-(M+1)])

n0 = 500 # number of sampled locations
n_new = 50 # to make predictions
set.seed(231)
yind <- sample(1:(M^2), size=n0+n_new)

coords = griddf[yind,]
D = dist(coords, diag = T, upper = T) %>% as.matrix() %>% unname()
rownames(D) <- NULL
coords_new = tail(coords, n_new)

n_sim = 100
mu1 = rep(0, n0+n_new)
mu2 = rep(0, n0+n_new)

init = 7
end = 7
restore = F

type_dgm = "RelBias_deltafixed"
# 3 options are available => "conditional", "RelBias_deltafixed", "RelBias_sigma2fixed"


process_tprsMM <- function(mm, to_df=TRUE){
  # Reorder and drop intercept
  mm <- mm[, c(ncol(mm)-1, ncol(mm), 1:(ncol(mm)-3))]
  if (to_df) {
    mm <- as.data.frame(mm)
  } 
  colnames(mm) <- paste0("tprs", 1:ncol(mm))
  mm
}

options_model = 1 # null-space type: choose 1, 2, or 3



priors = list(beta = c(0, 1e6), # N(mean, var)
              sigma2y = c(2, 0.1), # IG(shape, rate)
              spikeslab = c(1e-4, 1), #(s_0, s_1)
              w = c(1, 1), #Beta(a_w, b_w)
              psi = c(2, 1) #IG(a_p, b_p)
)
SS_priortype = c("george"=1, "nmig"=2, "pMOM"=3)
prior_choice = SS_priortype["pMOM"] # choose prior structure

switch (names(prior_choice),
        'george' = {subfolder = paste0("SS", options_model)},
        'nmig' = {subfolder = paste0("SS_nmig", options_model)},
        'pMOM' = {
          subfolder = paste0("SS_mom", options_model)
          priorCoef=momprior(tau=0.348)
          priorDelta <- modelbbprior(alpha.p=1,beta.p=1)
          priorVar <- igprior(alpha = priors$sigma2y[1], lambda = priors$sigma2y[2])
          varpmom = integrate(function(x) x^2 * dmom(x,tau=priorCoef@priorPars[1]), -Inf, Inf)$value #standardized variance
          if (priorCoef@priorDistr=='pMOM') method <- 'auto' else method <- 'Laplace'
        }
)

setting2 = expand.grid(range1 = seq(0.05, 0.5, length.out=10),
                       range2 = seq(0.05, 0.5, length.out=10),
                       delta = 0.5, sigma_1 = 1, sigma_2 = 1, tau = 0.25)
setting2 = data.frame(sim=1:nrow(setting2), cov_func='exp',
                      smooth_param=NA, setting2)


find.param2 = function(object, names, thin){
  nchains = length(object)
  out = lapply(1:nchains, function(i) mcmc(object[[i]][,names]))
  return(as.mcmc.list(out, thin=thin))
}

cl <- makeCluster(2, outfile = "foo.txt")
registerDoParallel(cl)
mycomb <- function(x, ...) {
  mapply(rbind,x,...,SIMPLIFY=FALSE)
}

# iterations and burnin are small numbers for the purpose of illustration only
# results in the paper are based on much larger values
iterations = 100
burnin = 100
thinning = 1
priorCoef=momprior(tau=0.348)
priorDelta <- modelbbprior(alpha.p=1,beta.p=1)
priorVar <- igprior(alpha = priors$sigma2y[1], lambda = priors$sigma2y[2])


psplines = principal.spline.2D(coords/max(coords), 1:n0, ~lngcoords+latcoords)
B = as.data.frame(psplines$F_n[,-1])
tprsm <- c(3:9, seq(10, 50, by=5), seq(60, 90, by=10), seq(100, 250, by=25))
# CREATE TPRS BASIS MATRIX
tprsSC.grid <- smoothCon(s(lngcoords, latcoords, fx=TRUE, k=n0,
                           xt=list(max.knots=2e3, seed=1776)), data = coords[1:n0,])
tprsMM0 <- PredictMat(tprsSC.grid[[1]], data = coords[1:n0,])
tprsDF <- process_tprsMM(tprsMM0, to_df=F)

betaPMOM = betaOLS = betaKS = matrix(NA, nrow(setting2), n_sim)
betaSPPLUS_fx = betaSPPLUS = betaSPATIALTP = betaSRE = betaGSEM = betaSA = betaKS

for (j in setting2$sim) {
  print(j)
  
  DG = generate_data(setting2, j, type_dgm, coords, D, 0.15)
  Xsim = DG$X
  Ysim = DG$Y
  
  estsTPRS <- matrix(0, nrow=length(tprsm), ncol=n_sim)
  rownames(estsTPRS) <- tprsm
  varTPRS = bicyzTPRS = aicyzTPRS = s2TPRS = fvMSE_TPRS = estsTPRS

  for(i in 1:n_sim){
    X <- Xsim[i,1:n0]
    Y <- Ysim[i,1:n0]
    x_new = tail(Xsim[i,], n_new)
    y_new = tail(Ysim[i,], n_new)
    
    mod = lm(Y~X)
    betaOLS[j,i] = mod$coefficients[2]
    
    mae = vector(l=length(tprsm)); names(mae) = tprsm; mspe = mae
    y_pred = matrix(nrow = n_new, ncol = length(tprsm))
    colnames(y_pred) = tprsm
    se_pred = y_pred
    
    # KELLER-SZPIRO
    for (k in 1:length(tprsm)){
      Xtprs_semipar_temp <- data.frame(int=1, X=X, as.matrix(tprsDF[, paste0("tprs", 1:tprsm[k])])) # prima c'era tprsDF[yind,...] mentre ora ho messo i valori campionati
      mod = lm(Y~-1 + ., data = Xtprs_semipar_temp)
      assign(paste('modKS',i,k, sep='_'), mod)
      sm = summary(mod)
      estsTPRS[k, i] <- coef(mod)[ind]
    }
    # AIC for outcome model with no exposure
    for (k in 1:length(tprsm)){
      Ztprs_semipar_temp <- as.matrix(cbind(1, as.matrix(tprsDF[, paste0("tprs", 1:tprsm[k])]))) # prima c'era tprsDF[yind,...] mentre ora ho messo i valori campionati
      mod = lm(Y~-1+Ztprs_semipar_temp)
      aicyzTPRS[k, i] <- AIC(mod)
      bicyzTPRS[k, i] <- BIC(mod)
    }
    if(i%%5 == 0) cat("i=",i,"\n")
    
    # SS_MOM
    xa = model.matrix(~., cbind(X, B))[,-1]
    includevars <- rep(FALSE, ncol(xa))
    includevars[2:ncol(model.matrix(~., data.frame(X=X)))-1] = T
    mod_SS_pMOM = histmat_pMOM= vector("list", 1)
    nchains=1
    for (chain in 1:nchains) {
      tic()
      set.seed(NULL)
      mod_SS_pMOM[[chain]] = modelSelection(
        y=Y, x=xa, center=T, scale=T,
        includevars = includevars, niter=iterations+burnin,
        burnin = burnin, priorCoef=priorCoef,
        thinning = thinning, priorDelta=priorDelta, priorVar=priorVar, 
        method=method, initSearch='greedy', verbose=T)
      histmat_pMOM[[chain]] <- rnlp(
        msfit=mod_SS_pMOM[[chain]], priorCoef=priorCoef, priorVar=priorVar,
        niter = iterations/thinning, burnin = burnin, thinning = thinning)
      toc()
    }
    beta_mcmc = find.param2(histmat_pMOM, 2, thinning)
    s = summary(beta_mcmc)
    betaPMOM[j,i] = unname(s$statistics['Mean'])
    
    # SPATIALTP, SPATIAL+, GSEM
    mod<-gam(Y~X+s(lngcoords,latcoords,k=150,fx=F),data=coords[1:n0,],method='GCV.Cp')
    betaSPATIALTP[j,i] = mod$coefficients[2]
    
    model_fx = T
    mod_X <- gam(X~s(lngcoords,latcoords,k=150,fx=model_fx),data=coords[1:n0,],method='GCV.Cp')
    f_X_hat = mod_X$fitted.values
    r_X<-X-f_X_hat
    mod<-gam(Y~r_X+s(lngcoords,latcoords,k=150, fx=model_fx),data=coords[1:n0,],method='GCV.Cp')
    betaSPPLUS_fx[j,i] = mod$coefficients[2]
    
    model_fx = F
    mod_X <- gam(X~s(lngcoords,latcoords,k=150,fx=model_fx),data=coords[1:n0,],method='GCV.Cp')
    f_X_hat = mod_X$fitted.values
    r_X<-X-f_X_hat
    mod<-gam(Y~r_X+s(lngcoords,latcoords,k=150, fx=model_fx),data=coords[1:n0,],method='GCV.Cp')
    betaSPPLUS[j,i] = mod$coefficients[2]
    
    mod_Y <- gam(Y~s(lngcoords,latcoords,k=150,fx=model_fx),data=coords[1:n0,],method='GCV.Cp')
    f_Y_hat = mod_Y$fitted.values
    r_Y<-Y-f_Y_hat
    mod = lm(r_Y~r_X-1)
    betaGSEM[j,i] = mod$coefficients[1]
  }
  ind_ic_min <- apply(aicyzTPRS, 2, which.min)
  est <- estsTPRS[cbind(ind_ic_min, 1:ncol(estsTPRS))]
  betaKS[j,] = unname(est[1:n_sim])
  
  # SRE
  fe = foreach(i=1:n_sim,.combine = 'mycomb', .packages = c('spBayes')) %dopar% {
    X<-Xsim[i,1:n0]
    Y<- Ysim[i,1:n0]
    x_new = tail(Xsim[i,], n_new)
    y_new = tail(Ysim[i,], n_new)
    
    mod = spLM(Y ~ X, coords = as.matrix(coords[1:n0,]),
               starting = list(beta=c(0,0), sigma.sq=.1, tau.sq=.1, phi=1/.2, nu=.5),
               tuning = list(sigma.sq=.02, tau.sq=.05, phi=.01, nu=0),
               priors = list(sigma.sq.ig=c(2,1), tau.sq.ig=c(2,1),
                             beta.norm=list(c(0,0), diag(1e6,2)),
                             phi.unif=c(1/.5,1/.05)), verbose = F, n.report = 10,
               cov.model = 'exponential', n.samples = iterations+burnin)
    m <- spRecover(mod, start = burnin, thin = thinning, verbose = F)
    
    beta_results = colMeans(m$p.beta.recover.samples)[ind]
    if(i%%1 == 0) capture.output(cat("isim=",i,"\n"), file="foo.txt", append=T)
    list(beta_results=beta_results)
  }
  betaSRE[j,] = fe$beta_results
  
  # SA
  SA_results = runsim_spectral(Xsim, Ysim, iterations, burnin, thinning, beta.real, 
                               predict_ = F, ncores = 2, log_file = 'log_SAmodel.txt')
  betaSA[j,] = SA_results$pmeanbeta[,2]
  
  
}


stopCluster(cl)
closeAllConnections()

all.res = list(setting2=setting2, betaPMOM = betaPMOM, betaKS=betaKS,
               betaSPPLUS_fx = betaSPPLUS_fx, betaOLS = betaOLS,
               betaSPPLUS = betaSPPLUS, betaSPATIALTP = betaSPATIALTP,
               betaSRE = betaSRE, betaGSEM = betaGSEM, betaSA = betaSA)
saveRDS(all.res, paste0(outpath, '/ratio.models.OLS.RDS'))
















