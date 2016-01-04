
gaussianmix <- function(datas, G, iter, weights, eqlam, tol = 1e-3, sigma.control = 1e-3, key = "kernelGaussian", 
                        verbose = FALSE, cov.constr = TRUE, cal.mallowslike = TRUE, mc = 0.25)
{
  # @param datas instances in rows
  # @param sigma.control upscales too small sigma estimate
  # @param cov.constr set covariance matrix diagonal if TRUE, general matrix if FALSE
  # @param verbose print det of covariance if TRUE
  # @param cal.mallowslike calculate mallows loglike if TRUE
  
  mnorm <- function(x, cent, covar, nconst = NULL, use.logC = TRUE){
    # @param x matrix of observations with samples in rows and features in col
    # @param cent a vector of mean
    # @param covar a single variance parameter such that covariance matrix is covar*IndentityMatrix
    # @param nconst normalization constant
    # @return log of multi-variate gaussian density for each row in x with mean cent and covariance matrix covar*I
    
    if (is.null(nconst)) {
      nconst <- -log(2*pi*covar)*ncol(x)/2
    } else if (!use.logC) {
      nconst <- log(nconst)
    }
    
    nconst - apply(t(x)-cent, 2, crossprod)/(2*covar)
  }
  
  updateProb <- function(N, G, x, mu, sigma, alpha, cov.constr){
    if (cov.constr) {
      p <- lapply(1:G, function(i) mnorm(x, cent = mu[[i]], covar = sigma[[i]]))
      p <- t(do.call("rbind", p) + log(alpha))
    } else {
      p <- lapply(1:N, function(j){
        log(alpha) + sapply(1:G, function(i){
          dmvnorm(x[j, ], mean = mu[[i]], sigma = sigma[[i]], log = TRUE)
        })
      })
      p <- do.call("rbind", p)
    }
    logsum <- logsumexp(p, byrow = TRUE)
    exp(p - logsum)
  }
  
  updateAlpha <- function(p, N, weights){
    colSums(p * weights) / sum(weights)
  }
  
  updateMu <- function(i, x, p, weights){
    pcol <- p[ ,i] * weights
    colSums(pcol * x) / sum(pcol)
  }
  
  updateSigma <- function(i, x, p, weights, mu, dfs, cov.constr, sigma.control){
    # @param sigma.control used to control fairly small sigma estimate when restricted model is employed
    # return for non-equal sigma estimate
    
    pcol <- p[ ,i] * weights
    if (cov.constr) {
      s <- sum(apply(t(x) - mu[[i]], 2, crossprod) * pcol) / sum(pcol) / dfs
      max(s, sigma.control, na.rm = TRUE)
    } else {
      covar <- apply(t(x) - mu[[i]], 2, tcrossprod)
      dim(covar) <- c(dfs, dfs, nrow(x))
      apply(aperm(covar, c(3,1,2))*pcol, c(2,3), sum) / sum(pcol)
    }
  }
  
  updateSigma_Eqlam <- function(x, p, weights, mu, dfs, cov.constr, sigma.control){
    # return for equal sigma estimate
    
    s <- lapply(1:G, function(i) apply(t(x) - mu[[i]], 2, crossprod))
    s <- do.call("cbind", s)
    pw <- weights * p
    s <- sum(s * pw) / sum(pw) / dfs
    max(s, sigma.control, na.rm = TRUE)
  }
  
  gaussianlike <- function(x, N, G, mu, alpha, sigma, weights, cov.constr, nconsts = NULL, use.logC = TRUE){
    if (cov.constr) {
      ff <- lapply(1:G, function(i) mnorm(x, cent = mu[[i]], covar = sigma[[i]], nconst = nconsts[[i]], use.logC = use.logC))
      ff <- logsumexp(do.call("rbind", ff) + log(alpha), bycol = TRUE)
    } else {
      ff <- sapply(1:N, function(j){
        probs <- sapply(1:G, function(i){
          dmvnorm(x[j, ], mean = mu[[i]], sigma = sigma[[i]], log = TRUE)
        })
        logsumexp(log(alpha) + probs)
      })
    }
    sum(weights * ff)
  }
  
  if (eqlam && !cov.constr) stop("equal covariance matrices that are not diagonal - not implemented yet")
  
  x <- KendallInfo(datas)
  N <- nrow(x)
  dfs <- ncol(x)
  
  # init
  p <- matrix(runif(N*G, min = 0, max = 1), N, G); p <- p/rowSums(p)
  alpha <- rep(1/G, G)
  mu <- lapply(1:G, function(i) runif(dfs, min = -1, max = 1))
  sigma <- rep(list(runif(1, min = 0.1, max = 10)), G)
  if (!cov.constr) sigma <- lapply(sigma, function(s) s*diag(ncol(x)))
  if (verbose) cat("iter = ", 0, ", det(SIGMA) = ", if(cov.constr) unlist(sigma)^dfs else sapply(sigma, det), "\n")
  like <- gaussianlike(x = x, N = N, G = G, mu = mu, alpha = alpha, sigma = sigma, weights = weights, cov.constr = cov.constr)
  
  likelihood <- 0*(1:iter)
  i <- 1
  
  while(i <= iter){
    # E-step
    p <- updateProb(N = N, G = G, x = x, mu = mu, sigma = sigma, alpha = alpha, cov.constr = cov.constr)
    
    # M-step
    alpha <- updateAlpha(p = p, N = N, weights = weights)
    mu <- lapply(1:G, updateMu, x = x, p = p, weights = weights)
    if (eqlam) {
      sigma <- updateSigma_Eqlam(x = x, p = p, weights = weights, mu = mu, dfs = dfs, cov.constr = cov.constr, sigma.control = sigma.control)
      sigma <- rep(list(sigma), G)
    } else {
      sigma <- lapply(1:G, updateSigma, x = x, p = p, weights = weights, mu = mu, dfs = dfs, cov.constr = cov.constr, sigma.control = sigma.control)
    }
    if (verbose) cat("iter = ", i, ", det(SIGMA) = ", if(cov.constr) unlist(sigma)^dfs else sapply(sigma, det), "\n")
    
    like <- gaussianlike(x = x, N = N, G = G, mu = mu, alpha = alpha, sigma = sigma, weights = weights, cov.constr = cov.constr)
    
    likelihood[i] <- like
    
    if (i > 2 && is.finite(likelihood[i]) && is.finite(likelihood[i-1]) && abs(likelihood[i] - likelihood[i-1]) < tol) {
      # message("Algorithm converged")
      likelihood[i:iter] <- likelihood[i]
      i <- iter
    }
    
    i <- i + 1
  }
  
  all.dists.data <- distL2(r = x, centers = do.call("rbind", mu))
  # lambda = 1/sigma only defined in restricted case; NOTE d(r,R) = ||r.info - R.info||^2/2 so that lambda = 1/sigma
  lambda <- if (cov.constr) 1/(2*unlist(sigma))/mc else rep(NA, length(sigma))
  out <- FormatOut(R = mu, p = alpha, lambda = lambda, z = p, datas = datas, likelihood = likelihood, 
                   all.dists.data = all.dists.data, weights = weights, key = key)
  
  if (cov.constr && cal.mallowslike) {
    abils <- ncol(datas)
    perm <- do.call("rbind", permn(abils))
    perm.info <- KendallInfo(perm)
    # search for true normalization constant such that it forms a distribution over permutation group
    nconsts <- lapply(1:G, function(i){
      -logsumexp(mnorm(x = perm.info, cent = mu[[i]], covar = sigma[[i]], nconst = 0, use.logC = TRUE))
    })
    # here is Mallows likelihood while keeping the estimate of all parameters
    out$min.like.mallows <- gaussianlike(x = x, N = N, G = G, mu = mu, alpha = alpha, sigma = sigma, weights = weights, cov.constr = cov.constr, nconsts = nconsts, use.logC = TRUE)
  }
  
  out
}
