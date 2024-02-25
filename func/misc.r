principal.spline.2D = function(coords, subset, formula, type = "radial", .sill = NULL, .range = NULL, .nu = 0.5, ...) {
  # coords: dataframe with spatial coordinates
  # subset: numeric vector with row indices indicating which rows of data are used to create the basis matrix, and which are used to create the prediction matrix
  # formula: a 'formula' object tha specifies how the trend fields are constructed
  # type: the potential function to use -> must be "radial" (d^2*log(d) where d is Euclidean distance)
  #                                        or "matern" (covariance function with parameters .sill, .range, .nu)
  # .sill: partial sill
  # .range: effective range parameter
  # .nu: smoothness parameter
  # ...: to pass other variables that are used in the formula
  
  N = length(subset)
  if(all(all.vars(formula) %in% names(coords))) {
    data = coords } else {
      args = list(...)
      data = cbind(coords, args)
      stopifnot("'formula' cannot find the variables. Check names of the arguments." = all.vars(formula) %in% names(data))
  }
  U.n = model.matrix(formula, data[subset,])
  U.n = unname(U.n); rownames(U.n) = NULL
  q = NCOL(U.n)
  coord.noti = coords[subset,] # coordinates of sampled sites

  potf = function(M, type) { # potential function
    M2 = switch(type,
                "radial" = M^2 * log(M),
                "matern" = cov.spatial(M, cov.pars = c(.sill, .range), kappa = .nu)
    )
    M2[is.nan(M2)] = 0
    return(M2)
  }

  dist2m = function(coo) { # equivalent to fields::rdist(coo)
    ans = unname(as.matrix(dist(coo)))
    rownames(ans) = NULL
    return(ans)
  }
  
  D.n = dist2m(coord.noti)
  S.n = potf(D.n, type) # NxN
  
  Ka = cbind(S.n, U.n) # Nx(N+q)
  Kb = cbind(t(U.n), 0*diag(q)) # q x (N+q)
  K = rbind(Ka, Kb) #(N+q)x(N+q)
  invK = solve(K)
  B.n = invK[1:N, 1:N] # bending energy matrix
  
  # decomposition bending energy matrix in space
  SVD = svd(B.n)
  FF = SVD$u # matrix whose columns contain the left singular vectors of B.n
  Ls = SVD$d # vector containing the singular values, sorted decreasingly
  FF = FF[, order(Ls)] # ascending order of left vectors
  Ls = sort(Ls) # ascending order of singular values
  Ls = Ls[-(1:q)]
  Fr1 = FF[, -(1:q)] # drop first q left vectors, corresponding to null singular values
  F.n = cbind(U.n, Fr1)
  
  if(N < nrow(data)) {
    N.p = nrow(data) - N
    coord.prev = coords[-subset,] # coordinates of new sites
    U.p = model.matrix(formula, data[-subset,])
    U.p = unname(U.p); rownames(U.p) = NULL
    
    D.p = dist2m(coord.prev)
    
    Dfull = dist2m(rbind(coord.noti, coord.prev))
    D.pn = Dfull[1:N, (N+1):ncol(Dfull)] # distances among each element in coord.prev to each in coord.noti
    D.pn = t(D.pn)
    
    S.p = potf(D.p, type)
    S.pn = potf(D.pn, type)
    A = invK[1:N, (N+1):ncol(invK)]
    U.p = unname(as.matrix(U.p)); rownames(U.p) = NULL
    alpha = tcrossprod(B.n, S.pn) + tcrossprod(A, U.p)
    F.p = cbind(U.p, crossprod(alpha, Fr1)) # matrix of bases for predictions in new sites
    
    res = list(F_n = F.n, F_p = F.p, q = q)
  } else {res = list(F_n = F.n, q = q)}
  
  return(res)
}

principal.spline.1D = function(X, subset, formula, type = "abs3", ...) {
  # X: dataframe with one column (e.g. representing time points)
  # subset: numeric vector with row indices indicating which rows of data are used to create the basis matrix, and which are used to create the prediction matrix
  # formula: a 'formula' object tha specifies how the trend fields are constructed
  # type: the potential function to use -> must be "abs3" (|d|^3 where d is pairwise difference of elements in X)
  # ...: to pass other variables that are used in the formula
  
  stopifnot("X is not a dataframe" = class(X) == "data.frame")
  N = length(subset)
  if(all(all.vars(formula) %in% names(X))) {
    data = X } else {
      args = list(...)
      data = cbind(X, args)
      stopifnot("'formula' cannot find the variables. Check names of the arguments." = all.vars(formula) %in% names(data))
    }
  U.n = model.matrix(formula, data[subset,, drop=F])
  U.n = unname(U.n); rownames(U.n) = NULL
  q = NCOL(U.n)
  X.noti = X[subset,] # coordinates of sampled sites
  
  potf = function(M, type) { # potential function
    M2 = switch(type,
                "abs3" = abs(M)^3
    )
    M2[is.nan(M2)] = 0 
    return(M2)
  }
  
  D.n = outer(X.noti, X.noti, "-") # time differences
  S.n = potf(D.n, type) # NxN
  
  Ka = cbind(S.n, U.n) # Nx(N+q)
  Kb = cbind(t(U.n), 0*diag(q)) # q x (N+q)
  K = rbind(Ka, Kb) #(N+q)x(N+q)
  invK = solve(K)
  B.n = invK[1:N, 1:N] # bending energy matrix
  
  # decomposition bending energy matrix in space
  SVD = svd(B.n)
  FF = SVD$u # matrix whose columns contain the left singular vectors of B.n
  Ls = SVD$d # vector containing the singular values, sorted decreasingly
  FF = FF[, order(Ls)] # ascending order of left vectors
  Ls = sort(Ls) # ascending order of singular values
  Ls = Ls[-(1:q)]
  Fr1 = FF[, -(1:q)] # drop first q left vectors, corresponding to null singular values
  F.n = cbind(U.n, Fr1)
  
  if(N < nrow(data)) {
    N.p = nrow(data) - N
    X.prev = X[-subset,] # times to predict
    U.p = model.matrix(formula, data[-subset,, drop=F])
    U.p = unname(U.p); rownames(U.p) = NULL
    
    D.p = outer(X.prev, X.prev, "-")
    
    X.pn = c(X.noti, X.prev)
    Dfull = outer(X.pn, X.pn, "-")
    D.pn = Dfull[1:N, (N+1):ncol(Dfull)] # time differences among each element in X.prev to each in X.noti
    D.pn = t(D.pn)
    
    S.p = potf(D.p, type)
    S.pn = potf(D.pn, type)
    A = invK[1:N, (N+1):ncol(invK)]
    U.p = unname(as.matrix(U.p)); rownames(U.p) = NULL
    alpha = tcrossprod(B.n, S.pn) + tcrossprod(A, U.p)
    F.p = cbind(U.p, crossprod(alpha, Fr1)) # matrix of bases for predictions in new sites
    
    res = list(F_n = F.n, F_p = F.p, q = q)
  } else {res = list(F_n = F.n, q = q)}
  
  return(res)
}


# Delta_OLS in the paper
bias_ols = function(sigma2, cxz, z, sigma1, L2i1) {
  z1 = cbind(1,z)
  H = solve(crossprod(z1)) %*% t(z1)
  bias = cxz * sqrt(sigma2/sigma1) * H %*% L2i1 %*% z
  return(bias[2])
}

# Delta_GLS in the paper
bias_gls = function(sigma2, cxz, R2, sigmay, z, sigma1, L2i1) {
  z1 = cbind(1,z)
  Swz = sigma2 * (1 - cxz^2) * R2
  Syz = Swz + (sigmay + beta.real[ind]^2 * 0.01) * diag(n0)
  Syz_inv = solve(Syz)
  H = solve(crossprod(z1, Syz_inv) %*% z1, crossprod(z1, Syz_inv))
  bias = cxz * sqrt(sigma2/sigma1) * H %*% L2i1 %*% z
  return(bias[2])
}


g.ls = function(formula, data, invOmega) {
  # modified version of File src/library/stats/R/lm.R
  # invOmega: inverse of var-cov matrix of disturbance term
  cl = mf = match.call()
  m <- match(c("formula", "data", "invOmega"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  x <- model.matrix(mt, mf)
  
  P = chol(invOmega)
  y = P %*% y
  x = P %*% x
  
  z <- lm.fit(x, y)
  class(z) <- c(if(NCOL(y)>1) "mlm", "lm")
  z$na.action <- attr(mf, "na.action")
  z$xlevels <- .getXlevels(mt, mf)
  z$call <- cl
  z$terms <- mt
  z$model <- mf
  return(z)
}



# function taken from this repo: https://github.com/emilelatour/lamisc.git
count_duplicates <- function(data, ...) {

  # unquoted names for NSE calls, need quoted names separately for messages +
  # warnings
  uq_names <- as.list(substitute(list(...)))[-1L]
  df_name <- deparse(substitute(data))

  if (length(uq_names) == 0) {
    # if called on an entire data.frame with no specified variable names
    var_names <- names(data)
    nms <- rlang::syms(var_names)
    #message("No variable names specified - using all columns.\n")
  } else {
    # nms <- rlang::quos(X2:X3) %>%
    #   rlang::quos_auto_name()
    # var_names <- names(nms)
    nms <- rlang::enquos(...)
  }

  df <- data %>%
    dplyr::select(!!! nms)

  x <- do.call('paste', c(df, sep = '\r'))
  ox <- order(x)
  rl <- rle(x[ox])
  cbind(df[ox[cumsum(rl$lengths)], , drop = FALSE], dupe_count = rl$lengths)

}

round10 = function(x) {
  count = 0
  s = sign(x)
  x2 = abs(x)
  while (x2>10) {
    x2 = x2/10
    count = count+1
  }
  if (x2>5) {
    ans = 10^(count+1)
  } else {
    ans = 10^count
  }
  return(ans*s)
}



### functions for KS model ####
getEsts <- function(X, y, w=rep(1, length(y)), ind=2, extra_params=0){
  n <- length(y)
  p <- ncol(X)
  Xw <- X * w
  XtX <- crossprod(X, Xw)
  lmX <- solve(XtX, crossprod(Xw, y))
  res2 <- as.vector((y - X %*% lmX)^2)
  varX <- solve(XtX, t(X))
  varX <-  varX  %*% (w^2*res2*t(varX))
  ll <- 0.5 * (sum(log(w)) - n * (log(2 * pi) + 1 - log(n) + log(sum(w*res2))))
  aic <- -2*ll + 2*(p + extra_params + 1) # +1 for variance parameter
  bic <- -2*ll + log(n)*(p + extra_params + 1) # +1 for variance parameter
  if (length(ind)>1){
    out <- list(coef=lmX[ind], var=varX[ind, ind], ll=ll, aic=aic[1], bic=bic[1],
                s2 = sum(res2)/(n-p), fv_MSE = mean(res2))
  } else {
    out <- c(coef=lmX[ind], var=varX[ind, ind], ll=ll, aic=aic[1], bic=bic[1],
             s2 = sum(res2)/(n-p), fv_MSE = mean(res2))
  }
  out
}

process_tprsMM <- function(mm, to_df=TRUE){
  # Reorder and drop intercept
  mm <- mm[, c(ncol(mm)-1, ncol(mm), 1:(ncol(mm)-3))]
  if (to_df) {
    mm <- as.data.frame(mm)
  } 
  colnames(mm) <- paste0("tprs", 1:ncol(mm))
  mm
}



### functions for models in Dupont et al. (2022) ####
# available at: https://doi.org/10.1111/biom.13656
# NOTE: some have been slightly modified to organize the output differently
fit_models<-function(X,Y,lngcoords, latcoords, sim_data, k_sp=300,model_fx=TRUE, pred_data){
  n0 = length(Y)
  beta_hat<-rep(NA,4)
  names(beta_hat)<-c("Spatial","RSR","gSEM","Spatial+")
  fv_MSE=s2mat=cover=aic=bic=r2adj=dev_prop=crossval<-beta_hat
  pred = list()
  #Spatial model
  mod<-gam(Y~X+s(lngcoords,latcoords,k=k_sp,fx=model_fx),data=sim_data,method=method)
  beta_hat[1]<-mod$coefficients[2]
  fv_MSE[1]<-mean((mod$fitted.values-Y)^2) #MSE of fitted values
  s2mat[1] = mod$sig2
  aic[1] = AIC(mod); bic[1] = BIC(mod)
  ci <- beta_hat[1] + sqrt(mod$Vp[2,2]) %o% c(-1, 1) * qnorm(0.975)
  cover[1] = ci[1] < beta.real[ind] & ci[2] > beta.real[ind]
  pred_mod = predict.gam(mod, pred_data, se = T)
  y_pred = data.frame(value = pred_mod$fit, se = pred_mod$se.fit)
  mae = mean(abs(y_pred$value - pred_data$Y))
  mspe = mean((y_pred$value - pred_data$Y)^2)
  pred[[1]] = list(summaries = c(mae=mae, mspe=mspe), values = y_pred)
  ssum = summary(mod)
  r2adj[1] = ssum$r.sq
  dev_prop[1] = ssum$dev.expl
  #RSR model
  mod_list<-gam(Y~X+s(lngcoords,latcoords,k=k_sp,fx=model_fx),data=sim_data,fit=FALSE)
  B_sp<-mod_list$X[,-(1:2)]
  X1 = cbind(1,X)
  P <- X1 %*% solve(crossprod(X1)) %*% t(X1)
  B_sp_tilde<-(diag(n0)-P)%*%B_sp
  mod_list$X[,-(1:2)]<-B_sp_tilde
  mod<-gam(G=mod_list,method=method)
  beta_hat[2]<-mod$coefficients[2]
  fv_MSE[2]<-mean((mod$fitted.values-Y)^2)
  s2mat[2]=mod$sig2
  aic[2] = AIC(mod); bic[2] = BIC(mod)
  ci <- beta_hat[2] + sqrt(mod$Vp[2,2]) %o% c(-1, 1) * qnorm(0.975)
  cover[2] = ci[1] < beta.real[ind] & ci[2] > beta.real[ind]
  pred_mod = predict.gam(mod, pred_data, se = T)
  y_pred = data.frame(value = pred_mod$fit, se = pred_mod$se.fit)
  mae = mean(abs(y_pred$value - pred_data$Y))
  mspe = mean((y_pred$value - pred_data$Y)^2)
  pred[[2]] = list(summaries = c(mae=mae, mspe=mspe), values = y_pred)
  #gSEM
  #Identify spatial patterns
  mod_X <- gam(X~s(lngcoords,latcoords,k=k_sp,fx=model_fx),data=sim_data,method=method)
  f_X_hat = mod_X$fitted.values
  r_X<-X-f_X_hat
  r_X_new = pred_data$X - predict.gam(mod_X, pred_data)
  mod_Y <- gam(Y~s(lngcoords,latcoords,k=k_sp,fx=model_fx),data=sim_data,method=method)
  f_Y_hat = mod_Y$fitted.values
  r_Y<-Y-f_Y_hat
  f_Y_new = predict.gam(mod_Y, pred_data)
  #Fit linear model to residuals
  mod<-lm(r_Y~-1+r_X)
  beta_hat[3]<-mod$coefficients[1]
  fv_MSE[3]<-mean((mod$fitted.values+f_Y_hat-Y)^2)
  estimated.coefs = (n0-mod_X$df.residual) + (n0-mod_Y$df.residual) + 1
  s2mat[3]=crossprod(mod$residuals)/(n0 - estimated.coefs)
  aic[3] = AIC(mod) - 2*2 + 2*(estimated.coefs+1); bic[3] = BIC(mod) - 2*log(n0) + (estimated.coefs+1)*log(n0)
  var_est = s2mat[3] * solve(crossprod(r_X))
  var_est = drop(var_est)
  ci <- beta_hat[3] + sqrt(var_est) %o% c(-1, 1) * qnorm(0.975)
  cover[3] = ci[1] < beta.real[ind] & ci[2] > beta.real[ind]
  pred_mod = predict(mod, data.frame(r_X = r_X_new))
  y_pred = data.frame(value = pred_mod+f_Y_new, se = NA)
  mae = mean(abs(y_pred$value - pred_data$Y))
  mspe = mean((y_pred$value - pred_data$Y)^2)
  pred[[3]] = list(summaries = c(mae=mae, mspe=mspe), values = y_pred)
  ssum = summary(mod)
  r2adj[3] = ssum$adj.r.squared
  dev_prop[3] = deviance(mod)/deviance(glm(r_Y~1))
  #Spatial plus
  mod<-gam(Y~r_X+s(lngcoords,latcoords,k=k_sp, fx=model_fx),data=sim_data,method=method)
  beta_hat[4]<-mod$coefficients[2]
  fv_MSE[4]<-mean((mod$fitted.values-Y)^2)
  s2mat[4]=mod$sig2
  aic[4] = AIC(mod) + 2*(n0-mod_X$df.residual); bic[4] = BIC(mod) + (n0-mod_X$df.residual)*log(n0)
  ci <- beta_hat[4] + sqrt(mod$Vp[2,2]) %o% c(-1, 1) * qnorm(0.975)
  cover[4] = ci[1] < beta.real[ind] & ci[2] > beta.real[ind]
  pred_mod = predict.gam(mod, data.frame(pred_data, r_X = r_X_new), se = T)
  y_pred = data.frame(value = pred_mod$fit, se = pred_mod$se.fit)
  mae = mean(abs(y_pred$value - pred_data$Y))
  mspe = mean((y_pred$value - pred_data$Y)^2)
  pred[[4]] = list(summaries = c(mae=mae, mspe=mspe), values = y_pred)
  ssum = summary(mod)
  r2adj[4] = ssum$r.sq
  dev_prop[4] = ssum$dev.expl
  
  list(beta_hat=beta_hat,fv_MSE=fv_MSE, s2mat = s2mat, cover = cover, prediction = pred, aic = aic, bic = bic)
}





### functions for SA model ####
# from repo: https://github.com/yawenguan/spatial_confounding.git
# NOTE: some have been slightly modified to organize the output differently
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
                            iters=10000,burn=1000,update=1000,thin=1,Xgrid = NULL,block = F, init=NULL, MH=0.25){
  ti = proc.time()
  Y = matrix(Y,nc=1)
  lbx = NCOL(X)
  if(lbx==1) X = matrix(X,nc=1)
  
  # initial param for MCMC
  tau    <- ifelse(is.null(init), sum(resid(mod_ols)^2)/(n0-lbx), init$sigma2y)  # measurement error variance
  betaX  <- ifelse(is.null(init), coef(mod_ols), init$beta)
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
      VVV  <- solve(1/tau*x.inv.x + diag(0.001,lbx))
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
        # MH     <- 0.25 # sd of proposal density
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
          # MH   <- 0.25
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
      cat(iter*thin, "= ")
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

