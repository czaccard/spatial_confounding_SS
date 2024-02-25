par_f = function(i, SS_priortype, Ysim, includevars, iterations, burnin, thinning, priors, w_indifference, numBlocks, true_value, complete){

  cat("###### Simulation number:", i, "######\n")
  y = Ysim[i, 1:n0]
  y_new = tail(Ysim[i,], n_new)
  mod_ols = lm(y ~ x)
  
  w = sum50 = sum75 = which50 = which75 = NA
  
  set.seed(NULL)
  if (SS_priortype!=3) {
    initial = list()
    initial$beta = solve(crossprod(xb), crossprod(xb,y))
    initial$sigma2y = crossprod(y - xb %*% initial$beta) / (nrow(xb) - ncol(xb))
    initial$sigma2y = drop(initial$sigma2y)
    
    ris = spsl(y, xb, B2, initial, numBlocks, priors, SS_priortype, iterations, burnin, thinning, w_indifference = w_indifference, verbose = F)
    
    betas = cbind(ris$intercept_history, ris$Beta.draws)
    pmeanbeta = colMeans(betas)
    qtl = apply(betas, 2, quantile, probs = c(0.025, 0.975))
    qtl = unname(qtl)
    cilow = qtl[1,]
    cihigh = qtl[2,]
    pmeanbasis = colMeans(ris$xi.draws)
    priorvarxi = ris$Is_slab_matrix*ris$psi.draws + (1-ris$Is_slab_matrix)*priors$spikeslab[1]*ris$psi.draws
    pmeanpsi2 = colMeans(ris$psi.draws)
    pmeansigma2 = colMeans(ris$sigma2y.draws)
    w = colMeans(ris$w.draws)
    pip = colMeans(ris$Is_slab_matrix)
    
    hist_mat = cbind(betas, ris$xi.draws, ris$sigma2y.draws)
  } else {
    ms <- modelSelection(y=y, x=design[,-1], center=T, scale=T, includevars = includevars[-1], niter=iterations+burnin, burnin = burnin, priorCoef=priorCoef,
                         thinning=thinning, priorDelta=priorDelta, priorVar=priorVar, method=method, initSearch='greedy', verbose=F)
    hist_mat <- rnlp(msfit=ms, priorCoef=priorCoef, priorVar=priorVar, niter = iterations/thinning, burnin = burnin, thinning = thinning)
    postsummary <- data.frame(selection=c(1,ms$postMode), margpp=c(0,ms$margpp), mean=colMeans(hist_mat)[1:p], t(apply(hist_mat[,1:p],2,quantile,probs=c(.025,.975))))
    
    pmeanbeta = postsummary$mean[1:lbx]
    cilow = postsummary$X2.5.[1:lbx]
    cihigh = postsummary$X97.5.[1:lbx]
    pmeanbasis = postsummary$mean[-(1:lbx)]
    pmeansigma2 = mean(hist_mat[,"phi"])
    pmeanpsi2 = NA
    pip = postsummary$margpp[-(1:lbx)]

    priorvarxi = ms$postSample[,-(1:(lbx-1))] * varpmom * sd(y)^2 / varB + (1 - ms$postSample[,-(1:(lbx-1))]) * 1e-10 # intercept not included in postSample
  }
  truth = rep(NA, lbx); truth[1:length(true_value)] = true_value
  ci.coverage = (truth >= cilow & truth <= cihigh)
  # Effective sample size
  ess = hist_mat[,2] %>% coda::mcmc() %>% coda::effectiveSize()
  
  #Prediction
  nc = ncol(hist_mat)
  y_pred = apply(hist_mat, 1, function(k) {
    rnorm(length(y_new), drop(design_new %*% k[-nc]), sqrt(k[nc]))
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

  # WAIC & DIC
  ll = apply(hist_mat, 1, function(r) {
    dnorm(y, drop(design %*% r[-nc]), sqrt(r[nc]), log = T)
  })
  waic = unlist(LaplacesDemon::WAIC(ll))
  Dev <- -2*colSums(ll)
  dic <- c(DIC=mean(Dev) + var(Dev)/2, Dbar=mean(Dev), pV=var(Dev)/2)

  d_pmean = c(NA, NA)
  d_bias = NA
  if (complete) {
    d_pmean = c(0,0)
    # Posterior Mean vs. OLS estimate
    yhat = unname(predict(mod_ols))
    switch (as.character(options_model),
            '1' = {
              s0 = BTB_tBtemp2B / err_var
              dup_mat = count_duplicates(data.frame(priorvarxi))
              d_mcmc = apply(dup_mat %>% dplyr::select(!dupe_count), 1, function(v) {
                invVxi = diag(v^-1)
                s = s0 + invVxi
                invs = solve(s)
                t = temp1B %*% invs
                d = t %*% t(B2) %*% (yhat-y)/err_var
                return(d[2])
              })
              d_pmean[1] = weighted.mean(d_mcmc, dup_mat$dupe_count)
              d_mcmc = apply(cbind(hist_mat[,nc], priorvarxi), 1, function(v) {
                invVxi = diag(v[-1]^-1)
                s = BTB_tBtemp2B/v[1] + invVxi
                invs = solve(s)
                t = temp1B %*% invs
                d = t %*% t(B2) %*% (yhat-y)/v[1]
                return(d[2])
              })
              d_pmean[2] = mean(d_mcmc)
            },
            '2' = {
              d = temp1C %*% (yhat-y)
              d_pmean = rep(d[2], 2)
            },
            '3' = {d_pmean = rep(0,2)}
    )

    # Bias Analysis (OLS without adjustment vs. OLS adjusted for confounding using principal splines)
    switch (as.character(options_model),
            '1' = {
              dup_mat = count_duplicates(data.frame(if (SS_priortype!=3) ris$Is_slab_matrix else ms$postSample[,-(1:(lbx-1))]))
              d_mcmc = apply(dup_mat %>% dplyr::select(!dupe_count), 1, function(r) {
                r = as.logical(r)
                bb = B2[, r, drop=F]
                if (NCOL(bb)!=0) {
                  s = crossprod(bb) - t(bb) %*% temp2 %*% bb
                  invs = solve(s)
                  t = temp1 %*% bb %*% invs
                  temp4 = const * t %*% t(bb)
                  p1 = temp4 %*% temp2 %*% temp3
                  p2 = temp4 %*% temp3
                  return(p1[2] - p2[2])
                } else {return(0)}
              })
              d_bias = weighted.mean(d_mcmc, dup_mat$dupe_count)
            },
            '2' = {
              d_bias = temp423[2]
            },
            '3' = {
              d_bias = 0
            }
    )
  }
  
  out = list(pmeanbeta = pmeanbeta, cilow = cilow, cihigh = cihigh, ci.coverage = ci.coverage, pmeanbasis = pmeanbasis, pip = pip,
             pmeansigma2 = pmeansigma2, pmeanpsi2 = pmeanpsi2, d_pmean = d_pmean, d_bias = d_bias, ess = ess, waic = waic, dic = dic,
             w = w, prediction = prediction)
  
  return(out)
}

runParallel_sim = function(SS_priortype, Xsim, Ysim, options_model, iterations, burnin, thinning, priors, w_indifference = T,
                           numBlocks, true_value, complete = T, ncores, log_file = 'mylog.txt'){
  # w_indifference: the parameter w is set to 0.5 if TRUE (indifference between spike and slab), or a Beta prior is assumed if FALSE
  # numBlocks: to implement the block-update strategy for the coefficients, this should be greater than 1
  require(parallel)
  on.exit(stopCluster(cl))
  if (file.exists(log_file)) writeLines("",log_file) #Clear file if it exists
  
  stopifnot(options_model %in% 1:3)
  switch(as.character(options_model),
         '1' = {
           includevars = rep(F, n0+1); includevars[1:2] = T
           reg_formula = ~x+lngcoords+latcoords+B
           p_formula = ~lngcoords+latcoords
         },
         '2' = {
           includevars = rep(F, n0); includevars[1:4] = T
           reg_formula = ~x+lngcoords+latcoords+B
           p_formula = ~x+lngcoords+latcoords
         },
         '3' = {
           includevars = rep(F, n0); includevars[1:2] = T
           reg_formula = ~x+B
           p_formula = ~x
         }
  )
  stopifnot(all(apply(Xsim, 2, function(x) length(unique(x)) == 1))) # just to make sure that all rows are equal (the exposure is fixed in our simulations)
  n_sim = nrow(Xsim)
  p = length(includevars); lbx = sum(includevars)
  ngood = floor(iterations/thinning)
  
  call = match.call()
  call$SS_priortype=SS_priortype; call$n_sim=n_sim; call$iterations=iterations; call$burnin=burnin; call$thinning=thinning; call$priors = priors
  call$options_model = options_model

  x = Xsim[1, 1:n0]
  x_new = tail(Xsim[1,], n_new)
  
  psplines = principal.spline.2D(coords, 1:n0, p_formula, x = Xsim[1,])
  dim_null = psplines$q
  B = psplines$F_n[,-(1:dim_null)]
  B_new = psplines$F_p[,-(1:dim_null)]
  xa = cbind(1, x)
  reg_formula = formula(paste(str_replace(as.character(reg_formula), "B", paste0("B.", 1:ncol(B), collapse = "+")), collapse = " "))
  design_full = model.matrix(reg_formula, data.frame(x=c(x, x_new), coords, B=rbind(B, B_new)))
  design_full = unname(design_full)
  stopifnot(ncol(design_full) == length(includevars))
  design = design_full[1:n0,]
  design_new = tail(design_full, n_new)
  xb = design[, includevars]
  B2 = design[, !includevars]
  varB = matrix(1, nrow = ngood) %*% apply(B2, 2, var)
  BTB = crossprod(B2)
  temp1 = solve(crossprod(xa)) %*% t(xa)
  temp2 = xa %*% temp1
  temp3 = DG$Lw[1:n0,1:n0] %*% solve(DG$Lz[1:n0,1:n0]) %*% x
  C = unname(as.matrix(coords[1:n0,]))
  invcmat = C %*% solve(crossprod(C, diag(n0) - temp2) %*% C) %*% t(C)
  temp1C = temp1 %*% invcmat
  const = DG$delta * sqrt(DG$sigma2w/DG$sigma2z)
  temp4 = const * temp1C
  temp423 = temp4 %*% (temp2 - diag(n0)) %*% temp3
  expected_betaOLS = true_value[2] + const * solve(crossprod(xa)) %*% t(xa) %*% temp3
  temp1B = temp1 %*% B2
  tBtemp2B = t(B2) %*% temp2 %*% B2
  BTB_tBtemp2B = BTB - tBtemp2B
  err_var = DG$sigma2y # true nugget
  
  cl <- makePSOCKcluster(ncores, outfile = log_file)
  setDefaultCluster(cl)
  clusterCall(cl, function() {library(dplyr)})
  if(SS_priortype !=3) {
    clusterExport(cl, c('spsl', 'n0', 'n_new', 'count_duplicates'))
  } else {
    clusterExport(cl, c('modelSelection', 'rnlp', 'priorCoef', 'priorVar', 'priorDelta', 'method', 'varpmom',
                        'n0', 'n_new', 'count_duplicates'))}
  #clusterExport(cl, 'x', 'xb', envir = environment())
  clusterExport(cl, ls(envir = environment()), envir = environment())

  simR <- parLapply(cl, 1:n_sim, par_f, SS_priortype=SS_priortype, Ysim=Ysim, includevars=includevars, iterations=iterations, burnin=burnin,
                    thinning=thinning, priors=priors, w_indifference=w_indifference, numBlocks=numBlocks, true_value=true_value, complete=complete)
  
  out = list()
  out$pmeanbeta = t(sapply(1:n_sim, function(i) simR[[i]]$pmeanbeta)); colnames(out$pmeanbeta) = paste0("beta", 1:lbx-1)
  out$cilow = t(sapply(1:n_sim, function(i) simR[[i]]$cilow)); colnames(out$cilow) = paste0("beta", 1:lbx-1)
  out$cihigh = t(sapply(1:n_sim, function(i) simR[[i]]$cihigh)); colnames(out$cihigh) = paste0("beta", 1:lbx-1)
  out$ci.coverage = t(sapply(1:n_sim, function(i) simR[[i]]$ci.coverage)); colnames(out$ci.coverage) = paste0("beta", 1:lbx-1)
  out$pmeanbasis = t(sapply(1:n_sim, function(i) simR[[i]]$pmeanbasis))
  out$pmeanpsi2 = t(sapply(1:n_sim, function(i) simR[[i]]$pmeanpsi2))
  out$pip = t(sapply(1:n_sim, function(i) simR[[i]]$pip))
  out$pmeansigma2 = sapply(1:n_sim, function(i) simR[[i]]$pmeansigma2)
  out$w = sapply(1:n_sim, function(i) simR[[i]]$w)
  out$prediction = lapply(1:n_sim, function(i) simR[[i]]$prediction)
  out$d_pmean = t(sapply(1:n_sim, function(i) simR[[i]]$d_pmean)); colnames(out$d_pmean) = c("true_nugget", "post_nugget")
  out$d_bias = sapply(1:n_sim, function(i) simR[[i]]$d_bias)
  out$ess = sapply(1:n_sim, function(i) simR[[i]]$ess)
  out$waic = t(sapply(1:n_sim, function(i) simR[[i]]$waic)); colnames(out$waic) = c("WAIC", "lppd", "pWAIC", "pWAIC1")
  out$dic = t(sapply(1:n_sim, function(i) simR[[i]]$dic)); colnames(out$dic) = c("DIC", "Dbar", "pV")
  out$call = call
  
  cat("Parallel execution finished\n")
  return(out)
}

# Assume that the log file is named 'mylog.txt'
# FOR WINDOWS USERS ONLY (optional):
# to see the output of runParallel_sim(), open PowerShell in your working dir.
# (example: in Win 11 open file explorer, go to the working dir., right click and choose open in Terminal)
# Run these commands in PowerShell:
# powershell Clear-Content mylog.txt
# powershell Get-Content mylog.txt -wait
# Once R has finished executing the calculations, type Ctrl+C in PowerShell to interrupt the listening command

# FOR LINUX USERS ONLY (optional):
# to see the output of runParallel_sim(), open terminal in your working dir. and run the command:
# tail -f mylog.txt
# Once R has finished executing the calculations, type Ctrl+C in terminal



run_OLS_fast = function(vec) {
  ols = vector(l = length(vec))
  jj=1
  for (j in vec) {
    DG = generate_data(setting, j, type_dgm, coords, D, 0.15)
    Xsim = DG$X
    Ysim = DG$Y
    
    B = vector(l=n_sim)
    
    for(i in 1:n_sim) {
      X<-Xsim[i,1:n0]
      Y<- Ysim[i,1:n0]
      XD = cbind(1, X)
      mod<-lm(Y~-1 + XD)
      B[i]<-mod$coefficients[ind]
    }
    ols[jj] = mean(B)
    jj = jj+1
    if(jj%%10 == 0)cat(".")
  }
  cat("\n")
  return(ols)
}


