# in this file:
# g:      unmeasured confounder
# X:      exposure
# Y:      outcome


generate_data = function(data, j, type, coords, D, relbias = NULL) {
  # data: file excel
  # j: configuration number
  # type: data generating mechanism
  # D: distance matrix based on coords
  # relbias: amount of relative bias when this should be fixed
  options(warn=1)
  on.exit(options(warn=0))
  
  e = list(); alph = NULL
  e$call_DG_function = match.call()
  e$call_DG_function$j = j; e$call_DG_function$type = type; e$call_DG_function$relbias = relbias
  
  N = nrow(coords)
  sigmay = data$tau[j] #variance of measurement error (nugget)
  rho1 = data$range1[j] #range
  sigma1 = data$sigma_1[j] #partial sill
  cxz = data$delta[j]
  sigma2 = data$sigma_2[j]
  rho2 = data$range2[j]
  corfuntype = data$cov_func[j]
  smoothness.par = data$smooth_param
  if (corfuntype == "exp") {
    corfuntype = "exponential"
    smoothness.par = 0.5
  }
  alph = smoothness.par + 1
  if (alph>2) { warning("Changing alpha to 2"); alph = 2 }
  stopifnot(any(corfuntype == c("matern", "exponential", "gaussian", "spherical", "powered.exponential")))
  if (any(corfuntype == c("matern", "powered.exponential"))) {
    R1 = cov.spatial(D, cov.model = corfuntype, cov.pars = c(1, rho1), kappa = smoothness.par)
    R2 = cov.spatial(D, cov.model = corfuntype, cov.pars = c(1, rho2), kappa = smoothness.par)
  } else {
    R1 = cov.spatial(D, cov.model = corfuntype, cov.pars = c(1, rho1))
    R2 = cov.spatial(D, cov.model = corfuntype, cov.pars = c(1, rho2))
  }
  L1 = t(chol(R1)) # lower triangular
  L2 = t(chol(R2)) # lower triangular
  
  
  switch (type,
          "conditional" = {
            set.seed(5324)
            z = rmvnorm(1, mu1, sigma1 * R1, method = "chol") %>% t()
            
            # g
            mg = cxz*sqrt(sigma2/sigma1) * L2 %*% solve(L1)
            Sg = sigma2*(1-cxz^2)*R2
            set.seed(334)
            Gsim <- rmvnorm(n_sim, mu2 + mg %*% (z - mu1), Sg, method = "chol")
            
            # X
            set.seed(5462)
            X = z + rnorm(N, sd = 0.1)
            XD = cbind(1, X) # design matrix
            Xsim = matrix(rep(X, n_sim), nrow = n_sim, byrow = T)
            
            # Y
            Ysim = Xsim*0
            Xb = XD%*%beta.real
            set.seed(642)
            for(i in 1:n_sim){
              Ysim[i,] =  Xb + Gsim[i,] + rnorm(N, sd = sqrt(sigmay))
            }
          },
          "RelBias_deltafixed" = {
            if(is.null(relbias)) stop("Relbias is required")
            if(sign(cxz) != sign(relbias)) {
              warning("Changing sign to 'relbias'.")
              relbias = -1*relbias
            }
            
            set.seed(5324)
            z = rmvnorm(1, mu1, sigma1 * R1, method = "chol") %>% as.vector()
            L2i1 = L2 %*% solve(L1)
            
            sigma2 = uniroot(function(x) {
              bias_ols(x, cxz, z, sigma1, L2i1)/beta.real[ind] - relbias
            }, interval = c(0, 1e10))$root
            
            cat("partial sill of unmeasured confounder:", sigma2, "\n")
            
            # g
            mg = cxz*sqrt(sigma2/sigma1) * L2i1
            Sg = sigma2*(1-cxz^2)*R2
            
            set.seed(334)
            Gsim <- rmvnorm(n_sim, mu2 + mg %*% (z - mu1), Sg, method = "chol")
            
            # X
            set.seed(5462)
            X = z + rnorm(N, sd = 0.1)
            XD = cbind(1, X) # design matrix
            Xsim = matrix(rep(X, n_sim), nrow = n_sim, byrow = T)
            
            # Y
            Ysim = matrix(nrow = n_sim, ncol = N)
            Xb = XD%*%beta.real
            set.seed(642)
            for(i in 1:n_sim){
              Ysim[i,] =  Xb + Gsim[i,] + rnorm(N, sd = sqrt(sigmay))
            }
            
          },
          "RelBias_sigma2fixed" = {
            if(is.null(relbias)) stop("Relbias is required")
            if(sign(cxz) != sign(relbias)) {
              warning("Changing sign to 'relbias'.")
              relbias = -1*relbias
            }
            
            set.seed(5324)
            z = rmvnorm(1, mu1, sigma1 * R1, method = "chol") %>% as.vector()
            L2i1 = L2 %*% solve(L1)
            
            cxz = uniroot(function(x) {
              bias_ols(sigma2, x, z, sigma1, L2i1)/beta.real[ind] - relbias
            }, interval = c(-1, 1))$root
            cat("correlation:", cxz, "\n")
            stopifnot(cxz<.8 && cxz>-.8)
            
            # g
            mg = cxz*sqrt(sigma2/sigma1) * L2i1
            Sg = sigma2*(1-cxz^2)*R2
            
            set.seed(334)
            Gsim <- rmvnorm(n_sim, mu2 + mg %*% (z - mu1), Sg, method = "chol")
            
            # X
            set.seed(5462)
            X = z + rnorm(N, sd = 0.1)
            XD = cbind(1, X) # design matrix
            Xsim = matrix(rep(X, n_sim), nrow = n_sim, byrow = T)
            
            # Y
            Ysim = matrix(nrow = n_sim, ncol = N)
            Xb = XD%*%beta.real
            set.seed(642)
            for(i in 1:n_sim){
              Ysim[i,] =  Xb + Gsim[i,] + rnorm(N, sd = sqrt(sigmay))
            }
          }
  )
  out = list(X = Xsim, Y = Ysim, alpha = alph, G = Gsim, e = e, sigma2y = sigmay)
  out2 = list(R_phiz=R1, R_phiw=R2, Lz=L1, Lw=L2, delta=cxz, sigma2z=sigma1, sigma2w=sigma2)
  out = c(out, out2)
  return(out)
}
