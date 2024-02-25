par_sp = function(i, Ysim, iterations, burnin, thinning, true_value){

  cat("###### Simulation number:", i, "######\n")
  y = Ysim[i, 1:n0]
  y_new = tail(Ysim[i,], n_new)
  mod_ols = lm(y ~ x)

  force2sym=function(s) { # needed because of some rare numerical fuzzy in the Prediction phase that prevents the variance matrix from being symmetric
    s[lower.tri(s)] = t(s)[lower.tri(s)]
    return(s)
  }
  
  vecc <- function(r,rho,nu,S,nn){
    # covparms for vecchia_Linv: see NOTE in this help page
    # ?GpGp::matern_isotropic
    
    Li  <- vecchia_Linv(c(r,rho,nu,1/r-1),
                        "matern_isotropic",
                        S,nn)
    return(Li)}
  # log(det(inverse_cov))
  vecc_logdet <- function(Li){
    2*sum(log(Li[,1]))
  }
  # t(y)%*%inverse_cov%*%y for vector y
  vecc_qf_vec <- function(Li,y,nn){
    sum(Linv_mult(Li,y,nn)^2)
  }
  # t(y)%*%inverse_cov%*%y for matrix y
  vecc_qf_mat <- function(Li,y,nn){
    Ly <- NULL
    for(j in 1:ncol(y)){
      Ly <- cbind(Ly,Linv_mult(Li,y[,j],nn))
    }
    out <- t(Ly)%*%Ly
    return(out)}

  # t(x)%*%inverse_cov%*%y
  vecc_mult <- function(x,Li,y,nn){
    if(is.vector(x)){
      Lx <- Linv_mult(Li,x,nn)
    }
    if(is.matrix(x)){
      Lx <- NULL
      for(j in 1:ncol(x)){
        Lx <- cbind(Lx,Linv_mult(Li,x[,j],nn))
      }
    }
    if(is.vector(y)){
      Ly <- Linv_mult(Li,y,nn)
    }
    if(is.matrix(y)){
      Ly <- NULL
      for(j in 1:ncol(y)){
        Ly <- cbind(Ly,Linv_mult(Li,y[,j],nn))
      }
    }
    out <- t(Lx)%*%Ly
    return(out)
  }
  SemiPcausal_vecc<- function(Y,X,d=NULL,L=100,nugget=TRUE,nn=30,zhat = NULL,
                    iters=10000,burn=1000,update=1000,thin=1,Xgrid = NULL,block = F){
    ti = proc.time()
    Y = matrix(Y,nc=1)
    lbx = NCOL(X)
    if(lbx==1) X = matrix(X,nc=1)
    
    # initial param for MCMC
    tau    <- sum(resid(mod_ols)^2)/(n0-2) # measurement error variance
    betaX  <- coef(mod_ols)
    betaL  <- rnorm(L)
    betaL = betaL - sum(betaL)/L
    rhores <- 0.3 # range
    nures  <- 0.5 # smoothness
    r      <- ifelse(nugget,0.5,1) # var/totalvar 
    tauprecb <- 1 # precision for ICAR prior for b_l
    Q      <- fields::rdist(1:L)==1 
    Q      <- diag(rowSums(Q))-Q # ICAR prior for b_l, coefficients of zhat
    # Q      <- diag(rowSums(Q))-0.99*Q # proper prior for b_l, coefficients of zhat
    
    # place holder
    keep.betaX <- matrix(0,iters,lbx)
    keep.betaL <- matrix(0,iters,L)
    keepers    <- matrix(0,iters,4)
    dev        <- rep(0,iters)
    colnames(keepers) <- c("tau","r","range","sig2b")
    accrates = matrix(0, iters, 2, dimnames = list(NULL, c("range", "r")))
    
    
    # precompute a few quantities
    n       <- length(Y)
    distmat <- fields::rdist(d)
    # compute confounder adjustment if not given
    # s        = Xgrid # assume X observed at a fine grid
    # alpha    = seq(0,30*pi,l=L)
    # if(is.null(zhat)) zhat = getzhat_gaussian_grid(rdist(d, s),m=alpha,rho=rhozhat,nu=nuzhat,X,rdist(s[1:2,])[2])
    nn_ordered <- find_ordered_nn(d,nn)
    Li       <- vecc(r,rhores,nures,d,nn_ordered) # SigmaRes = maternFun(distmat, c(r,rhores,nures,1/r-1))
    curdet   <- vecc_logdet(Li) # -determinant(SigmaRes)$mod
    
    
    for(iter in 1:iters){
      for(thinthin in 1:thin){
        accrho = accr = F
        # Update beta
        Yc       <- Y-zhat%*%betaL
        x.inv.x  <- vecc_qf_mat(Li,X,nn_ordered) # (t(cbind(1,X))%*%solve(SigmaRes,cbind(1,X)))
        x.inv.y  <- vecc_mult(X,Li,Yc,nn_ordered) # t(cbind(1,X))%*%solve(SigmaRes,Yc)
        VVV  <- solve(1/tau*x.inv.x + diag(0.0001,lbx))
        MMM  <- 1/tau*x.inv.y
        betaX<- VVV%*%MMM+t(chol(VVV))%*%rnorm(lbx)
        # Update betazhat
        Yc     <- Y-X%*%betaX
        zhat.inv.zhat  <- vecc_qf_mat(Li,zhat,nn_ordered) #(t(zhat)%*%solve(SigmaRes,zhat))
        zhat.inv.y  <- vecc_mult(zhat,Li,Yc,nn_ordered) # (t(zhat)%*%solve(SigmaRes,Yc))
        VVV    <- solve(1/tau*zhat.inv.zhat + tauprecb*Q)
        MMM    <- 1/tau*zhat.inv.y
        #M      <- VVV%*%MMM
        betaL  <- VVV%*%MMM+t(chol(VVV))%*%rnorm(L)
        betaL = betaL - sum(betaL)/L
        
        # # setting the last one to zero. Do you need this?
        # condM  <- M[-L] + VVV[1:(L-1),L]/VVV[L,L]*(0-M[L])
        # condV  <- VVV[1:(L-1),1:(L-1)] + VVV[1:(L-1),L]%*%t(VVV[1:(L-1),L])/VVV[L,L]
        # betaL[1:(L-1)] <- condM + t(chol(condV))%*%rnorm(L-1)

        # update total variance and prior precision
        res  <- Y-X%*%betaX - zhat%*%betaL
        SSE  <- vecc_qf_vec(Li,res,nn_ordered) #t(res)%*%solve(SigmaRes,res)
        SSB  <- t(betaL)%*%Q%*%betaL
        tau  <- 1/rgamma(1,n/2+2,SSE/2+1)
        tauprecb <- rgamma(1,L/2+2,SSB/2+1)
        # VVV  <- 1/(1/tau*tauprecb*sum(Q) + 0.01)
        # MMM  <- 1/tau*tauprecb*sum(Q%*%betaL)
        # mub  <- rnorm(1,VVV*MMM,sqrt(VVV)) # Forgot what this is for
        
        
        if(!block){
          # Update rhores, prior unif(0,1.5)
          MH     <- 0.25 # sd of proposal density
          canrho <- exp(rnorm(1,log(rhores),MH))
          canLi  <- vecc(r,canrho,nures,d,nn_ordered)# canSigmaRes <- maternFun(distmat, c(r,canrho,nures,1/r-1))
          candet <- vecc_logdet(canLi) # -determinant(canSigmaRes)$mod 
          canSSE <- vecc_qf_vec(canLi,res,nn_ordered) #t(res)%*%solve(canSigmaRes,res)
          
          canl <- 1/2*candet - 1/(2*tau)*canSSE + dunif(canrho, 0, 1.5, log = T)
          curl <- 1/2*curdet - 1/(2*tau)*SSE + dunif(rhores, 0, 1.5, log = T)
          R    <- canl-curl-log(rhores)+log(canrho)
          
          if(log(runif(1))<R){
            rhores   <- canrho
            Li       <- canLi # SigmaRes <- canSigmaRes
            curdet   <- candet
            SSE      <- canSSE# t(res)%*%solve(SigmaRes,res)
            accrho = T
          }
          
          # Update r, prior unif(0,1)
          if(nugget){
            MH   <- 0.25
            canr <- pnorm(rnorm(1,qnorm(r),MH))
            canLi  <- vecc(canr,rhores,nures,d,nn_ordered)# canSigmaRes <- maternFun(distmat, c(canr,canrho,nures,1/canr-1))
            candet <- vecc_logdet(canLi) # -determinant(canSigmaRes)$mod 
            canSSE <- vecc_qf_vec(canLi,res,nn_ordered) #t(res)%*%solve(canSigmaRes,res)
            
            canl <- 1/2*candet - 1/(2*tau)*canSSE
            curl <- 1/2*curdet - 1/(2*tau)*SSE
            R    <- canl-curl-dnorm(qnorm(r),log=TRUE)+dnorm(qnorm(canr),log=TRUE)
            
            if(log(runif(1))<R){
              r      <- canr
              Li     <- canLi 
              curdet <- candet
              SSE    <- canSSE
              accr = T
            }
          } 
        }else{
          #  need to check this block of code for block update
          MH     <- 0.25
          canrho <- exp(rnorm(1,log(rhores),MH))
          canr   <- pnorm(rnorm(1,qnorm(r),MH))
          canLi  <- vecc(canr,canrho,nures,d,nn_ordered)# canSigmaRes <- maternFun(distmat, c(r,canrho,nures,1/r-1))
          candet <- vecc_logdet(canLi) # -determinant(canSigmaRes)$mod 
          canSSE <- vecc_qf_vec(canLi,res,nn_ordered) #t(res)%*%solve(canSigmaRes,res)
          
          canl <- 1/2*candet - 1/(2*tau)*canSSE + dunif(canrho, 0, 1.5, log = T)
          curl <- 1/2*curdet - 1/(2*tau)*SSE + dunif(rhores, 0, 1.5, log = T)
          R    <- canl-curl-log(rhores)+log(canrho)-dnorm(qnorm(r),log=TRUE)+dnorm(qnorm(canr),log=TRUE)
          
          U = rnorm(2)
          
          
          if(log(runif(1))<R){
            rhores   <- canrho
            r        <- canr
            Li       <- canLi # SigmaRes <- canSigmaRes
            curdet   <- candet
            SSE      <- canSSE# t(res)%*%solve(SigmaRes,res)
          }
        }
      } # end thin
      
      if(iter%%update==0){
        capture.output(print(iter), file="logggg.txt", append=TRUE)
      }
      
      
      dev[iter]         <- -2*(-n/2*log(tau) + 1/2*curdet - 1/(2*tau)*SSE)
      keep.betaX[iter,] <- betaX
      keep.betaL[iter,] <- betaL
      keepers[iter,]    <- c(tau,r,rhores,1/tauprecb)
      accrates[iter,] = c(accrho, accr)
    } # end iter
    
    #parameter means
    mn       <- colMeans(keepers[(burn+1):iters,])
    betaX.mn <- colMeans(keep.betaX[(burn+1):iters,])
    betaL.mn <- colMeans(keep.betaL[(burn+1):iters,])
    accrates = accrates[(burn+1):iters,]
    res.mn <- Y-X%*%betaX.mn - zhat%*%betaL.mn
    tau.mn <- mn[1]
    Li.mn  <- vecc(mn[2],mn[3],nures,d,nn_ordered)# canSigmaRes <- maternFun(distmat, c(r,canrho,nures,1/r-1))
    det.mn <- vecc_logdet(Li.mn) # -determinant(canSigmaRes)$mod 
    SSE.mn <- vecc_qf_vec(Li.mn,res.mn,nn_ordered) #t(res)%*%solve(canSigmaRes,res)
    
    Dbar <- mean(dev[(burn+1):iters])
    Dhat <- -2*(-n/2*log(tau.mn ) + 1/2*det.mn - 1/(2*tau.mn)*SSE.mn)
    pD   <- Dbar-Dhat
    DIC  <- Dbar + pD
    DIC  <- list(DIC=DIC,Dbar=Dbar,pD=pD)
    time = proc.time() - ti
    out <- list(betaX=keep.betaX[(burn+1):iters,],betaL=keep.betaL[(burn+1):iters,],dev=dev,keepers=keepers[(burn+1):iters,],
                DIC = DIC,time=time, accrates = accrates)
    return(out)
  }

  alldic = NULL
  deltas = seq(1,40,by=5)
  for(del in deltas){
    zhat = getzhat_bspline_grid(coords, region2[,-3], L, X=region2$z, del=del,setup = setup) #calculate confounder adjustment zhat
    set.seed(NULL)
    fitsp_bs = SemiPcausal_vecc(y, cbind(1,x), d=as.matrix(coords[1:n0,]), L=ncol(zhat$zhat), update=1e7, nugget=F,
                                zhat=zhat$zhat[1:n0,],iters=(iterations+burnin)/thinning,
                                burn = burnin/thinning, thin = thinning, nn=20)
    alldic = c(alldic, fitsp_bs$DIC$DIC)
    
    assign(paste0("fitsp_bs_",del),fitsp_bs)
    assign(paste0("zhat_",del),zhat)
  }

  del = deltas[which.min(alldic)]
  ris = get(paste0("fitsp_bs_",del))
  zhat = get(paste0("zhat_",del))
  betas = ris$betaX
  pmeanbeta = colMeans(betas)
  qtl = apply(betas, 2, quantile, probs = c(0.025, 0.975))
  qtl = unname(qtl)
  cilow = qtl[1,]
  cihigh = qtl[2,]
  pmeanbasis = colMeans(ris$betaL)
  pmeansigma2 = mean(ris$keepers[,1])
  pmeanrho = mean(ris$keepers[,3])
  pmeansig2b = mean(ris$keepers[,4])
  hist_mat = cbind(betas, ris$betaL, ris$keepers[,3], ris$keepers[,1])

  truth = rep(NA, lbx); truth[1:length(true_value)] = true_value
  ci.coverage = (truth >= cilow & truth <= cihigh)
  # Effective sample size
  ess = hist_mat[,2] %>% coda::mcmc() %>% coda::effectiveSize()
  
  #Prediction
  nc = ncol(hist_mat)
  design = cbind(1, x, zhat$zhat[1:n0,])
  design_new = cbind(1, x_new, tail(zhat$zhat, n_new))
  y_pred = apply(hist_mat, 1, function(k) {
    gamma = k[nc]*exp(-dist2dist/k[nc-1])
    Sig = geoR::cov.spatial(distmat, cov.pars = c(k[nc], k[nc-1]))
    Sig_new = geoR::cov.spatial(distmat_new, cov.pars = c(k[nc], k[nc-1]))
    invS = solve(Sig)
    mvtnorm::rmvnorm(1, 
          drop(design_new %*% k[-c(nc-1,nc)] + gamma %*% invS %*% (y - design%*%k[-c(nc-1,nc)])), 
          force2sym(Sig_new - gamma %*% invS%*%t(gamma)), method = 'chol')
  })
  y_pred = apply(y_pred, 1, function(x) {
    r = data.frame(r=x)
    r = r %>% summarise(mean(r), quantile(r, 0.025), quantile(r, 0.975), crossprod(r))
    return(r)
  })
  y_pred = matrix(unlist(y_pred), ncol = 4, byrow = T, dimnames = list(NULL, c("mean", "q025", "q975", "ss"))) %>%
    as.data.frame() %>% mutate(var = ss/ngood - mean^2, sd = sqrt(var))
  mae = mean(abs(y_pred$mean - y_new))
  mspe = mean((y_pred$mean - y_new)^2)
  prediction = list(summaries = c(mae=mae, mspe=mspe), values = y_pred)

  # DIC
  dic <- c(DIC=ris$DIC$DIC, Dbar=ris$DIC$Dbar, pV=ris$DIC$pD)
  
  out = list(pmeanbeta = pmeanbeta, cilow = cilow, cihigh = cihigh, ci.coverage = ci.coverage, pmeanbasis = pmeanbasis, pmeanrho = pmeanrho,
             pmeansigma2 = pmeansigma2, pmeansig2b = pmeansig2b, ess = ess, dic = dic, prediction = prediction)
  
  return(out)
}

runsim_spectral = function(Xsim, Ysim, iterations, burnin, thinning, true_value, ncores, log_file){
  require(parallel)
  on.exit(stopCluster(cl))
  if (file.exists(log_file)) writeLines("",log_file) #Clear file if it exists
  
  setup_bspline_grid<-function(s1,s2,X,dx = 0.5){
    # dx = 0.5 or 1,should be divisible by del
    deltax   = diff(sort(unique(s2[,1]))[1:2])
    deltay   = diff(sort(unique(s2[,2]))[1:2])
    delta    = c(deltax,deltay)
    max.del.x = floor(pi/max(delta))
    del.x = seq(0,max.del.x,by = dx)
    
    distmat = fields::rdist(s1,s2)
    foo = lapply(del.x,function(xx) besselJ(distmat*xx, nu = 0)%*%X)
    foo = do.call(cbind,foo)
    return(list(del.x=del.x,foo = foo,delta = delta,dx = dx))
  }
  getzhat_bspline_grid<-function(s1,s2,L,X,del,setup=NULL){
    del.x = setup$del.x
    foo = setup$foo
    dx = setup$dx
    max.del.x = max(del.x)
    Kdelta = floor(max.del.x/del)
    
    deltax   = diff(sort(unique(s2[,1]))[1:2])
    deltay   = diff(sort(unique(s2[,2]))[1:2])
    delta.x    = c(deltax,deltay)
    
    # get bspline basis fn
    del.xx = c(del.x,seq(max.del.x+dx,(Kdelta+1)*del,by = dx))
    del.xx = c(rev(-seq(dx,3*del,by = dx)), del.xx)
    
    bsbasis = bs(del.xx,knots = seq(-3*del,(Kdelta+1)*del,by=del))
    # truncate to [0,wt]
    bsbasis = bsbasis[-c(which(del.xx<0),which(del.xx>max.del.x)),]
    del.xx  = del.xx[-c(which(del.xx<0),which(del.xx>max.del.x))]
    # keep only nonzero basis
    keep = colSums(bsbasis)!=0
    bsbasis = bsbasis[,keep]
    L = ncol(bsbasis)
    
    zhat = NULL
    for(l in 1:L){
      temp1 = 0
      for(xx in 1:length(del.xx)){
        temp1 = temp1 + 2*pi*del.xx[xx]*bsbasis[xx,l]*foo[,xx]
      }
      temp1 = temp1*diff(del.x)[1]*prod(delta.x)
      zhat = cbind(zhat,temp1)
    }
    zhat = zhat/(2*pi)^2
    return(list(zhat=zhat,del.xx=del.xx,bsbasis=bsbasis))
  }

  stopifnot(all(apply(Xsim, 2, function(x) length(unique(x)) == 1))) # just to make sure that all rows are equal (the exposure is fixed in our simulations)
  n_sim = nrow(Xsim)
  ngood = floor(iterations/thinning)
  lbx=2
  
  call = match.call()
  call$n_sim=n_sim; call$iterations=iterations; call$burnin=burnin; call$thinning=thinning

  region2 <- with(cbind(coords, x=Xsim[1,]), interp(x = lngcoords, y = latcoords, z = x, linear = FALSE, extrap = TRUE,
                                                xo=seq(0, max(griddf), length=M/2), 
                                                yo=seq(0, max(griddf), length=M/2)))
  region2 <- as.data.frame(interp2xyz(region2))
  setup = setup_bspline_grid(coords,region2[,-3],region2$z,dx = 1)

  x = Xsim[1, 1:n0]
  x_new = tail(Xsim[1,], n_new)
  
  distmat = fields::rdist(coords[1:n0,])
  distmat_new = fields::rdist(coords_new)
  dist2dist = fields::rdist(coords_new, coords[1:n0,])

  cl <- makePSOCKcluster(ncores, outfile = log_file)
  setDefaultCluster(cl)
  clusterCall(cl, function() {library(dplyr); library(splines); library(GpGp)})
  clusterExport(cl, c('n0', 'n_new', 'coords', 'coords_new', 'outpath'))
  clusterExport(cl, ls(envir = environment()), envir = environment())

  simR <- parLapply(cl, 1:n_sim, par_sp, Ysim=Ysim, iterations=iterations, burnin=burnin,
                    thinning=thinning, true_value=true_value)
  
  out = list()
  out$pmeanbeta = t(sapply(1:n_sim, function(i) simR[[i]]$pmeanbeta)); colnames(out$pmeanbeta) = paste0("beta", 1:lbx-1)
  out$cilow = t(sapply(1:n_sim, function(i) simR[[i]]$cilow)); colnames(out$cilow) = paste0("beta", 1:lbx-1)
  out$cihigh = t(sapply(1:n_sim, function(i) simR[[i]]$cihigh)); colnames(out$cihigh) = paste0("beta", 1:lbx-1)
  out$ci.coverage = t(sapply(1:n_sim, function(i) simR[[i]]$ci.coverage)); colnames(out$ci.coverage) = paste0("beta", 1:lbx-1)
  out$pmeanbasis = t(sapply(1:n_sim, function(i) simR[[i]]$pmeanbasis))
  out$pmeansig2b = sapply(1:n_sim, function(i) simR[[i]]$pmeansig2b)
  out$pmeansigma2 = sapply(1:n_sim, function(i) simR[[i]]$pmeansigma2)
  out$pmeanrho = sapply(1:n_sim, function(i) simR[[i]]$pmeanrho)
  out$prediction = lapply(1:n_sim, function(i) simR[[i]]$prediction)
  out$ess = sapply(1:n_sim, function(i) simR[[i]]$ess)
  out$dic = t(sapply(1:n_sim, function(i) simR[[i]]$dic)); colnames(out$dic) = c("DIC", "Dbar", "pV")
  out$call = call
  
  cat("Spectral Adjustment: Parallel execution finished\n")
  return(out)
}


# Assume that the log file is named 'mylog.txt'
# FOR WINDOWS USERS ONLY (optional):
# to see the output of runsim_spectral(), open PowerShell in your working dir.
# (example: in Win 11 open file explorer, go to the working dir., right click and choose open in Terminal)
# Run these commands in PowerShell:
# powershell Clear-Content mylog.txt
# powershell Get-Content mylog.txt -wait
# Once R has finished executing the calculations, type Ctrl+C in PowerShell to interrupt the listening command

# FOR LINUX USERS ONLY (optional):
# to see the output of runsim_spectral(), open terminal in your working dir. and run the command:
# tail -f mylog.txt
# Once R has finished executing the calculations, type Ctrl+C in terminal

