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

setting <- readxl::read_excel(paste0(inpath, "setting.xlsx"), sheet = "Foglio2")
#sourceCpp(paste0(mypath, "functions.cpp"))
source(paste0(inpath, "data_generate.r"), local = T)
source(paste0(inpath, "sim_functions.r"), local = T)
source(paste0(inpath, "misc.r"), local = T)

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


#### Spatial+ (DWA) ####
if (!dir.exists(paste0(outpath ,"DWA"))) dir.create(paste0(outpath ,"DWA"))


tprsm = kgrid <- c(3:9, seq(10, 50, by=5), seq(60, 90, by=10), seq(100, 250, by=25))
names_kgrid = paste0("tprs", tprsm)
model_names<-c("OLS","Spatial_fx","RSR_fx","gSEM_fx","Spatial+_fx","Spatial","RSR","gSEM","Spatial+")
method<-"GCV.Cp"


k_sp<-150
for (j in init:end) {
  DG = generate_data(setting, j, type_dgm, coords, D, 0.15)
  Xsim = DG$X
  Ysim = DG$Y
  
  beta_results<-matrix(0,length(model_names),n_sim)
  rownames(beta_results)<-model_names
  fv_MSE_results = coverage_results = s2_results<-beta_results
  appo = data.frame(matrix(NA, n_new, length(model_names)))
  colnames(appo) = model_names
  prediction = vector(l=n_sim) %>% as.list()
  prediction = lapply(1:n_sim, function(ii) prediction[[ii]] = list(summaries=matrix(NA, length(model_names), 2, dimnames=list(model_names, c('mae', 'mspe'))),
                                                                    values = appo, se = appo))
  
  for(i in 1:n_sim){
    X <- Xsim[i,1:n0]
    Y <- Ysim[i,1:n0]
    x_new = tail(Xsim[i,], n_new)
    y_new = tail(Ysim[i,], n_new)
    pred_data = data.frame(coords_new, X=x_new, Y=y_new)
    
    XD = cbind(1, X) # design matrix
    #Fit linear model
    mod<-lm(Y~X)
    beta_results[1,i]<-mod$coefficients[ind]
    fv_MSE_results[1,i]<-mean((mod$fitted.values-Y)^2)
    s2_results[1,i] = crossprod(mod$residuals)/(n0-ncol(XD))
    var_est = s2_results[1,i] * solve(crossprod(XD))[ind,ind]
    ci <- beta_results[1,i] + sqrt(var_est) %o% c(-1, 1) * qnorm(0.975)
    coverage_results[1,i] = ci[1] < beta.real[ind] & ci[2] > beta.real[ind]
    pred_mod = predict(mod, pred_data, se = T)
    y_pred = data.frame(value = pred_mod$fit, se = pred_mod$se.fit)
    prediction[[i]]$summaries[1,] = c(mae = mean(abs(y_pred$value - pred_data$Y)), mspe = mean((y_pred$value - pred_data$Y)^2))
    prediction[[i]]$values[,1] = y_pred[,1]
    prediction[[i]]$se[,1] = y_pred[,2]
    #Fit models without smoothing of spatial effects
    fit<-fit_models(X,Y, lngcoords, latcoords, sim_data = coords[1:n0,], k_sp=k_sp,model_fx=TRUE, pred_data = pred_data)
    beta_results[2:5,i]<-fit$beta_hat
    fv_MSE_results[2:5,i]<-fit$fv_MSE
    s2_results[2:5,i] = fit$s2mat
    coverage_results[2:5,i] = fit$cover
    prediction[[i]]$summaries[2:5,] = t(sapply(1:4, function(ii) fit$prediction[[ii]]$summaries))
    prediction[[i]]$values[,2:5] = sapply(1:4, function(ii) fit$prediction[[ii]]$values[,1])
    prediction[[i]]$se[,2:5] = sapply(1:4, function(ii) fit$prediction[[ii]]$values[,2])
    #Fit models with smoothing of spatial effects
    fit<-fit_models(X,Y, lngcoords, latcoords, sim_data = coords[1:n0,], k_sp=k_sp,model_fx=FALSE, pred_data = pred_data)
    beta_results[6:9,i]<-fit$beta_hat
    fv_MSE_results[6:9,i]<-fit$fv_MSE
    s2_results[6:9,i] = fit$s2mat
    coverage_results[6:9,i] = fit$cover
    prediction[[i]]$summaries[6:9,] = t(sapply(1:4, function(ii) fit$prediction[[ii]]$summaries))
    prediction[[i]]$values[,6:9] = sapply(1:4, function(ii) fit$prediction[[ii]]$values[,1])
    prediction[[i]]$se[,6:9] = sapply(1:4, function(ii) fit$prediction[[ii]]$values[,2])
    if(i%%10 == 0)cat("i=",i,"\n")
  }
  MSE_results = apply(beta_results,1, function(x) mean((x - beta.real[ind])^2))
  MSE_results = data.frame(mse = MSE_results)
  
  DWA_results = list(beta_results=t(beta_results), fv_MSE_results=t(fv_MSE_results), MSE_results=MSE_results, s2_results=t(s2_results),
                     prediction = prediction, coverage_results=t(coverage_results), beta.real=beta.real, input=data.frame(setting[j,]))
  saveRDS(DWA_results, file = paste0(outpath, "DWA/sim", j,"_", type_dgm,".RDS"))
  cat("Spatial+: finished simulations with setting # ", j, "\n")
}
rm(beta_results, fv_MSE_results, MSE_results, s2_results, coverage_results, var_est, ci)



#### Keller-Szpiro (KS) ####
if (!dir.exists(paste0(outpath ,"KS"))) dir.create(paste0(outpath ,"KS"))

tprsm <- c(3:9, seq(10, 50, by=5), seq(60, 90, by=10), seq(100, 250, by=25))
# CREATE TPRS BASIS MATRIX
tprsSC.grid <- smoothCon(s(lngcoords, latcoords, fx=TRUE, k=n0,
                      xt=list(max.knots=2e3, seed=1776)), data = coords[1:n0,])
tprsMM0 <- PredictMat(tprsSC.grid[[1]], data = coords[1:n0,])
tprsDF <- process_tprsMM(tprsMM0, to_df=F)
tprsMM0_new = PredictMat(tprsSC.grid[[1]], data = coords_new)
tprsDF_new <- process_tprsMM(tprsMM0_new, to_df=F)
# basiTPRS = cbind(1, tprsDF)

for (j in init:end) {
  DG = generate_data(setting, j, type_dgm, coords, D, 0.15)
  Xsim = DG$X
  Ysim = DG$Y
  
  estsTPRS <- matrix(0, nrow=length(tprsm), ncol=n_sim)
  rownames(estsTPRS) <- tprsm
  varTPRS = bicyzTPRS = aicyzTPRS = s2TPRS = fvMSE_TPRS = estsTPRS
  appo = list()
  
  for(i in 1:n_sim){
    X <- Xsim[i,1:n0]
    Y <- Ysim[i,1:n0]
    x_new = tail(Xsim[i,], n_new)
    y_new = tail(Ysim[i,], n_new)
    
    mae = vector(l=length(tprsm)); names(mae) = tprsm; mspe = mae
    y_pred = matrix(nrow = n_new, ncol = length(tprsm))
    colnames(y_pred) = tprsm
    se_pred = y_pred
    
    # TPRS-SemiPar
    for (k in 1:length(tprsm)){
      Xtprs_semipar_temp <- data.frame(int=1, X=X, as.matrix(tprsDF[, paste0("tprs", 1:tprsm[k])]))
      mod = lm(Y~-1 + ., data = Xtprs_semipar_temp)
      sm = summary(mod)
      estsTPRS[k, i] <- coef(mod)[ind]
      varTPRS[k, i] <- sm$coefficients[ind,2]^2
      s2TPRS[k, i] <- drop(crossprod(mod$residuals)/mod$df.residual)
      fvMSE_TPRS[k,i] = mean((Y - mod$fitted.values)^2) # MSE of fitted values
      pred_data = data.frame(1, x_new, tprsDF_new[, paste0("tprs", 1:tprsm[k])])
      colnames(pred_data) = colnames(Xtprs_semipar_temp)
      pred_mod = predict(mod, pred_data, se = T)
      y_pred[,k] = pred_mod$fit
      se_pred[,k] = pred_mod$se.fit
      mae[k] = mean(abs(pred_mod$fit - y_new))
      mspe[k] = mean((pred_mod$fit - y_new)^2)
    }
    appo[[i]] = list(summaries = data.frame(mae=mae, mspe=mspe), values = y_pred, se = se_pred)
    # AIC for outcome model with no exposure
    for (k in 1:length(tprsm)){
      Ztprs_semipar_temp <- as.matrix(cbind(1, as.matrix(tprsDF[, paste0("tprs", 1:tprsm[k])])))
      mod = lm(Y~-1+Ztprs_semipar_temp)
      aicyzTPRS[k, i] <- AIC(mod)
      bicyzTPRS[k, i] <- BIC(mod)
    }
    if(i%%5 == 0)cat("i=",i,"\n")
  }
  ind_ic_min <- apply(aicyzTPRS, 2, which.min)
  est <- estsTPRS[cbind(ind_ic_min, 1:ncol(estsTPRS))]
  estVar = varTPRS[cbind(ind_ic_min, 1:ncol(varTPRS))]
  ci <- est + sqrt(estVar) %o% c(-1, 1) * qnorm(0.975)
  
  KS_results = list()
  KS_results$beta_results <- est
  KS_results$s2_results <- s2TPRS[cbind(ind_ic_min, 1:ncol(s2TPRS))]
  KS_results$fv_MSE_results = fvMSE_TPRS[cbind(ind_ic_min, 1:ncol(fvMSE_TPRS))]
  KS_results$MSE_results <- mean((est - beta.real[ind])^2)
  KS_results$coverage_results = ci[,1] < beta.real[ind] & ci[,2] > beta.real[ind]
  KS_results$beta.real=beta.real
  KS_results$prediction = list(
    summaries = t(sapply(1:n_sim, function(ii) appo[[ii]]$summaries[ind_ic_min[ii],])),
    values = sapply(1:n_sim, function(ii) appo[[ii]]$values[,ind_ic_min[ii]]),
    se = sapply(1:n_sim, function(ii) appo[[ii]]$se[,ind_ic_min[ii]])
  )
  KS_results$input=data.frame(setting[j,])
  KS_results$notes = DG$e
  saveRDS(KS_results, file = paste0(mypath, "KS/sim", setting$sim[j],"_", type_dgm,".RDS"))
  cat("KS: finished simulations with setting # ", j, "\n")
}
rm(ind_ic_min, est, estVar, ci, appo, y_pred, se_pred, mae, mspe, mod, pred_mod)

#### Spike & Slab ####
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
if (!dir.exists(paste0(outpath , subfolder))) dir.create(paste0(outpath , subfolder))

iterations = 2000
burnin = 1000
thinning = 1

for (j in init:init) {
  DG = generate_data(setting, j, type_dgm, coords, D, 0.15)
  Xsim = DG$X
  Ysim = DG$Y
  
  SS_results = runParallel_sim(prior_choice, Xsim, Ysim, options_model, iterations, burnin, thinning, priors, T,
                       10, beta.real, T, ncores = 2)
  SS_results$input = data.frame(setting[j,])
  SS_results$MSE_results = mean((SS_results$pmeanbeta[,ind] - beta.real[ind])^2)
  SS_results$notes = DG$e
  savename = paste0(outpath, subfolder, "/sim", j,"_", type_dgm,".RDS")
  saveRDS(SS_results, file = savename)
  cat(subfolder, ": finished simulations with setting #", j, "\n")
}





#### Spectral Adjustment model ####
source(paste0(inpath, "sim_spectralAdj.r"), local = T)
subfolder = 'SpectralAdj'
iterations = 100
burnin = 50
thinning = 2

for (j in init:init) {
  DG = generate_data(setting, j, type_dgm, coords, D, 0.15)
  Xsim = DG$X
  Ysim = DG$Y
  
  SA_results = runsim_spectral(Xsim, Ysim, iterations, burnin, thinning,
                               beta.real, ncores = 2, log_file = 'log_SA.txt')
  SA_results$input = data.frame(setting[j,])
  SA_results$MSE_results = mean((SA_results$pmeanbeta[,ind] - beta.real[ind])^2)
  SA_results$notes = DG$e
  savename = paste0(outpath, subfolder, "/sim", j,"_", type_dgm,".RDS")
  saveRDS(SA_results, file = savename)
  cat(subfolder, ": finished simulations with setting #", j, "\n")
}





#### Spatial Random Effect ####
if (!dir.exists(paste0(outpath ,"SRE"))) dir.create(paste0(outpath ,"SRE"))

priors = list(beta = c(0, 1e6), # N(mean, var)
              sigma2y = c(2, 0.1), # IG(shape, rate)
              spikeslab = c(1e-4, 1), #(s_0, s_1)
              w = c(1, 1), #Beta(a_w, b_w)
              psi = c(2, 1) #IG(a_p, b_p)
)
iterations = 2000
burnin = 1000
thinning = 1

cl = makeCluster(3, outfile = 'log_SRE.txt')
registerDoParallel(cl)
mycomb <- function(x, ...) {
  mapply(rbind,x,...,SIMPLIFY=FALSE)
}

for (j in init:end) {
  capture.output(cat(j, "------"), file="log_SRE2.txt", append=T)
  DG = generate_data(setting, j, type_dgm, coords, D, 0.15)
  Xsim = DG$X
  Ysim = DG$Y
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
    m.p = spPredict(mod, start = burnin, thin = thinning, verbose = F,
                    pred.covars = cbind(1,x_new), pred.coords = as.matrix(coords_new))
    
    beta_results = colMeans(m$p.beta.recover.samples)[ind]
    thetas = cbind(m$p.theta.recover.samples, range=1/m$p.theta.recover.samples[,3]) # tau.sq=variance of nugget, sigma.sq=partial sill
    s2_results = colMeans(thetas[,-3])
    ci = apply(m$p.beta.recover.samples, 2, quantile, probs=c(2.5,97.5)/100)[,ind]
    coverage_results = ci[1] < beta.real[ind] & ci[2] > beta.real[ind]
    y_pred <- apply(m.p$p.y.predictive.samples, 1, mean)
    mae = mean(abs(y_pred - y_new))
    mspe = mean((y_pred - y_new)^2)
    if(i%%1 == 0) capture.output(cat("isim=",i,"\n"), file="foo.txt", append=T)
    list(beta_results=beta_results, coverage_results=coverage_results, s2_results= s2_results, mae=mae, mspe=mspe, y_pred=y_pred)
  }
  fe$MSE_results = mean((fe$beta_results - beta.real[ind])^2)
  fe$notes = DG$e
  
  savename = paste0(outpath, subfolder, "/sim", j,"_", type_dgm,".RDS")
  saveRDS(fe, file = savename)
  cat(subfolder, ": finished simulations with setting #", j, "\n")
}
stopCluster(cl)












