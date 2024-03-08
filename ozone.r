rm(list=ls())
outpath = "results/"
inpath = "func/"
plot_path = "plots/"
source(paste0(inpath, "misc.r"), local = T)

library(tidyverse)
library(spBayes)
library(mgcv)
library(SSconf)
library(mombf)
library(tictoc)
library(GpGp)
library(splines)
library(coda)


dat = readRDS("ozone_data.RDS")
coords = dat[,1:2]
stdevs = readRDS('stdevs.RDS')

#### OLS ####
f0 = reformulate(names(dat)[names(dat)!='y'], response = 'y')
f0 = update.formula(f0, ~.-1-lon-lat)
mod0 = lm(f0, dat)
summary(mod0)
coef(mod0)['x']
r = resid(mod0)
qqnorm(scale(r)[,1]); abline(0,1)

f = y~-1+x
mod_ols = lm(f, dat)
coef(mod_ols)['x']


#### Spatial: SRE ####
iterations = 2000
burnin = 1000
thinning = 1

mod0sp = spLM(f0, dat, as.matrix(coords),
              starting = list(beta=rnorm(7), sigma.sq=.1, tau.sq=.1, phi=1, nu=.5),
              tuning = list(sigma.sq=.04, tau.sq=.02, phi=.03, nu=0),
              priors = list(sigma.sq.ig=c(2,1), tau.sq.ig=c(2,1),
                            beta.norm=list(rep(0, 7), diag(1e6,7)),
                            phi.unif=c(3/max(dist(coords)), 3/min(dist(coords)))),
              cov.model = 'exponential', n.samples =iterations)
mod0sp = spRecover(mod0sp, start = burnin+1)
plot(mod0sp$p.beta.recover.samples[,'x'])

mod_sre = spLM(f, dat, as.matrix(coords),
               starting = list(beta=rnorm(1), sigma.sq=.1, tau.sq=.1, phi=1, nu=.5),
               tuning = list(sigma.sq=.02, tau.sq=.02, phi=.02, nu=0),
               priors = list(sigma.sq.ig=c(2,1), tau.sq.ig=c(2,1),
                             beta.norm=list(rep(0, 1), diag(1e6,1)),
                             phi.unif=c(3/max(dist(coords)), 3/min(dist(coords)))),
               cov.model = 'exponential', n.samples = iterations)
mod_sre = spRecover(mod_sre, start = burnin+1)
mean(mod_sre$p.beta.recover.samples[,'x'])



#### Spatial, Spatial+ and gSEM ####
n0 = nrow(dat)
k_sp = 200
fgam = update.formula(f, ~.+ s(lon, lat, k=k_sp, fx=T))
mod_spat_fx = gam(fgam, data = dat)
fgam = update.formula(f, ~.+ s(lon, lat, k=k_sp, fx=F))
mod_spat = gam(fgam, data = dat)

dat2 = select(dat, x, y, lon, lat)
mod_y <- gam(y ~ -1 + s(lon, lat,k=k_sp, fx=F), data=dat)
dat2$r_y <- dat2$y - mod_y$fitted.values
appo = select(dat, x)
for(cova in names(appo)){
  ff = reformulate('-1', cova)
  ff = update.formula(ff, ~.+ s(lon, lat,k=k_sp, fx=F))
  if(is.factor(get(cova, dat))) {
    fac2num = as.numeric(get(cova, dat))-1
    nCat = max(fac2num)
    fam = if(nCat==1) 'binomial' else {
      if(is.ordered(get(cova, dat))) quote(ocat(R=nCat+1)) else quote(multinom(nCat))
    }
  } else fam = 'gaussian'
  mod_c = gam(ff, data=dat2, family = eval(fam))
  i_c = which(names(dat2)==cova)
  dat2[,i_c] = dat2[,i_c] - mod_c$fitted.values
  names(dat2)[i_c] = paste0('r_', cova)
}

mod_gsem <- lm(r_y~.-1-y-lon-lat, data=dat2)
ff = reformulate(paste0('r_', names(appo)), response = 'y')
ff = update.formula(ff,~.-1+ s(lon, lat,k=k_sp, fx=F))
mod_plus <- gam(ff, data=dat2)




#### Keller-Szpiro ####
k_sp_max = round(n0/2,-2) # MAX basis size for spatial thin plate splines
tprsm <- c(5, 10, seq(0, k_sp_max, by=25)[-1])
appo = select(dat, x)
tprsSC.grid <- smoothCon(s(lon, lat, fx=T, k=k_sp_max,
                           xt=list(max.knots=2e3, seed=1776)), data = dat)
# tprsMM0 <- PredictMat(tprsSC.grid[[1]], data = dat)
tprsDF <- process_tprsMM(tprsSC.grid[[1]]$X, to_df=F)

saveAIC = saveBIC = NULL
for (k in 1:length(tprsm)){
  print(k/length(tprsm)*100)
  Xtprs_semipar_temp <- 
    data.frame(Y=dat$y, appo, as.matrix(tprsDF[, 1:(tprsm[k]-1)]))
  mod = lm(Y~ .-1, data = Xtprs_semipar_temp)
  assign(paste0("mod_KS_",k), mod)
  mod2 = lm(Y~ .-x-1, data = Xtprs_semipar_temp)
  saveAIC = c(saveAIC, AIC(mod2))
  saveBIC = c(saveBIC, BIC(mod2))
}

kstar = (1:length(tprsm))[which.min(saveAIC)]
mod_KS_AIC = get(paste0("mod_KS_",kstar))
kstar = (1:length(tprsm))[which.min(saveBIC)]
mod_KS_BIC = get(paste0("mod_KS_",kstar))
rm(list = paste0("mod_KS_",1:length(tprsm)))
# plot(variog(list(data=resid(mod_KS_BIC), coords=dat[,1:2]), 
#             uvec=seq(0,10,l=20)))
gc()


#### SS ####
appo = select(dat, x)
coords = dat %>% dplyr::select(lon, lat)
iterations = 1000
burnin = 1000
thinning = 1
nchains = 1
options_model = 1

priors = list(beta = c(0, 1), # N(mean, var)
              sigma2y = c(2, 0.1), # IG(shape, rate)
              spikeslab = c(1e-4, 1), #(s_0, s_1)
              w = c(1, 1), #Beta(a_w, b_w)
              psi = c(2, 1) #IG(a_p, b_p)
)
SS_priortype = c("george"=1, "nmig"=2, "pMOM"=3)

priorCoef=momprior(tau=0.348)
priorDelta <- modelbbprior(alpha.p=1,beta.p=1)
priorVar <- igprior(alpha = priors$sigma2y[1], lambda = priors$sigma2y[2])
varpmom = integrate(function(x) x^2 * dmom(x,tau=priorCoef@priorPars[1]), -Inf, Inf)$value #standardized variance
if (priorCoef@priorDistr=='pMOM') method <- 'auto' else method <- 'Laplace'

psplines = principal.spline.2D(coords/max(coords), 1:n0, ~lon+lat)
B = as.data.frame(psplines$F_n[,-1])
xa = model.matrix(~., cbind(appo, B))[,-1]
includevars <- rep(FALSE, ncol(xa))
includevars[2:ncol(model.matrix(~., appo))-1] = T
mod_SS_pMOM = histmat_pMOM= vector("list", nchains)
for (chain in 1:nchains) {
  tic()
  set.seed(NULL)
  mod_SS_pMOM[[chain]] = modelSelection(
    y=dat$y, x=xa, center=T, scale=T,
    includevars = includevars, niter=iterations+burnin,
    burnin = burnin, priorCoef=priorCoef,
    thinning = thinning, priorDelta=priorDelta, priorVar=priorVar, 
    method=method, initSearch='greedy', verbose=T)
  histmat_pMOM[[chain]] <- rnlp(
    msfit=mod_SS_pMOM[[chain]], priorCoef=priorCoef, priorVar=priorVar,
    niter = iterations/thinning, burnin = burnin, thinning = thinning)
  toc()
}
# plot(variog(list(data=dat$y-cbind(1,xa)%*%colMeans(histmat_pMOM[[1]][,1:(ncol(xa)+1)]), coords=dat[,1:2]), 
#             uvec=seq(0,10,l=20)))
mod_SS_FV = vector("list", nchains)
for (chain in 1:nchains) {
  tic()
  set.seed(NULL)
  mod_SS_FV[[chain]] = spsl(f, dat, B, NULL, 10, priors,
                            SS_priortype["george"], iterations, burnin, 
                            thinning, verbose = T, restore_work = F)
  toc()
}

mod_SS_NMIG = vector("list", nchains)
for (chain in 1:nchains) {
  tic()
  set.seed(NULL)
  mod_SS_NMIG[[chain]] = spsl(f, dat, B, NULL, 10, priors,
                              SS_priortype["nmig"], iterations, burnin, 
                              thinning, verbose = T, restore_work = F)
  toc()
}


#### Spectral Adjustment ####
iterations = 100
burnin = 100
thinning = 1
setup = setup_bspline_grid(as.data.frame(coords),as.data.frame(coords),dat$x, dx = 1)


initial = list()
alldic = NULL
deltas = seq(1,40,by=5)
for(del in deltas){
  print(paste(del, "-"))
  zhat = getzhat_bspline_grid(as.data.frame(coords),as.data.frame(coords), 
                              L, X=dat$x, del=del, setup = setup)
  initial$beta = rnorm(ncol(appo))
  initial$sigma2y = 1/rgamma(1,2)
  fitsp_bs = vector("list", nchains)
  for (chain in 1:nchains) {
    tic()
    set.seed(NULL)
    fitsp_bs[[chain]] = SemiPcausal_vecc(
      dat$y, cbind(1,as.matrix(appo)),
      d=as.matrix(coords), L=ncol(zhat$zhat), update=100, nugget=F,
      zhat=zhat$zhat,iters=(iterations+burnin)/thinning, init = initial,
      burn = burnin/thinning, thin = thinning, nn=20, MH = 1.5)
    toc()
  }
  
  alldic = c(alldic, mean(sapply(1:nchains, function(t) fitsp_bs[[t]]$DIC$DIC)))
  
  assign(paste0("fitsp_bs_",del),fitsp_bs)
  assign(paste0("zhat_",del),zhat)
}
delstar = deltas[which.min(alldic)]
print(delstar)
mod_SA=get(paste0("fitsp_bs_",delstar))
rm(list = paste0("fitsp_bs_",deltas))
sapply(1:nchains, function(t) colMeans(mod_SA[[t]]$accrates))
gc()



util_gam = function(mod, j) {
  if (is.character(j)) j = which(names(coef(mod))==j)
  ci_ = coef(mod)[j] + 
    sqrt(mod$Vp[j,j]) %o% c(-1, 1) * qnorm(0.975)
  return(c(coef(mod)[j], ci_))
}
find.param = function(object, name, thin){
  nchains = length(object)
  out = lapply(1:nchains, function(i) mcmc(get(name, object[[i]])))
  return(as.mcmc.list(out, thin=thin))
}
find.param2 = function(object, names, thin){
  nchains = length(object)
  out = lapply(1:nchains, function(i) mcmc(object[[i]][,names]))
  return(as.mcmc.list(out, thin=thin))
}


model_names = c("OLS", "Spatial+_fx", "SpatialTP", "gSEM", "Spatial+", "KS",
                "SS_fv", "SS_nmig", "SS_mom", "SA", "SRE")
table_x = data.frame(Model = model_names,
                     Covariate = "x",
                     Estimate = 0, LowCI = 0, UppCI = 0)
ci_ = confint(mod_ols)
table_x[1,3:5] = c(coef(mod_ols)['x'], ci_['x',])

table_x[2,3:5] = util_gam(mod_spat_fx, 'x')

table_x[3,3:5] = util_gam(mod_spat, 'x')

ci_ = confint(mod_gsem)
table_x[4,3:5] = c(coef(mod_gsem)['r_x'], ci_['r_x',])

table_x[5,3:5] = util_gam(mod_plus, 'r_x')

ci_ = confint(mod_KS_AIC)
table_x[6,3:5] = c(coef(mod_KS_AIC)['x'], ci_['x',])

beta_mcmc = find.param(mod_SS_FV, "Beta.draws", thinning)
s = summary(beta_mcmc[,'x'])
table_x[7,3:5] = c(s$statistics[1], s$quantiles[1], s$quantiles[5])

beta_mcmc = find.param(mod_SS_NMIG, "Beta.draws", thinning)
s = summary(beta_mcmc[,'x'])
table_x[8,3:5] = c(s$statistics[1], s$quantiles[1], s$quantiles[5])

beta_mcmc = find.param2(histmat_pMOM, 'x', thinning)
s = summary(beta_mcmc)
table_x[9,3:5] = c(s$statistics[1], s$quantiles[1], s$quantiles[5])

beta_mcmc = find.param(mod_SA, "betaX", thinning)
s = summary(beta_mcmc)
id_x = which(names(appo)=='x') +1
table_x[10,3:5] = c(s$statistics[id_x,1], s$quantiles[id_x,1], s$quantiles[id_x,5])

s = summary(mod_sre$p.beta.recover.samples[,'x'])
table_x[11,3:5] = c(s$statistics[1], s$quantiles[1], s$quantiles[5])


ci0x = confint(mod0)['x',]
table_x$low0=ci0x[1]
table_x$high0=ci0x[2]
table_x[,-(1:2)] = table_x[,-(1:2)]*stdevs['y']/stdevs['x']
model_names2 = c("OLS","SRE", "Spatial+_fx", "SpatialTP", "gSEM", "Spatial+", "KS",
                 "SA", "SS_fv", "SS_nmig", "SS_mom", "SGLASSO")
model_face = rep('plain', length(model_names2)-1)
model_face[9:11] = 'bold'


#### THIS IS FIGURE 3 ####
gg1rd = ggplot(table_x, aes(Model, Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin=LowCI, ymax=UppCI),
                position=position_dodge(0.05), width = 0.2, linewidth = .7) +
  geom_hline(aes(yintercept=0), linetype="dotted") +
  scale_x_discrete(limits=model_names2[-12]) +
  ylab('Estimated Effect (%)') +
  theme_light() +
  theme(strip.text = element_text(size=12, color="black"),
        panel.grid = element_blank(),
        axis.text.x = element_text(face = model_face, angle=90))
gg1rd
ggsave(paste0(plot_path, "application_ozone.pdf"), gg1rd,
       width = 7.5, height = 4, device = pdf)
