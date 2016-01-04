#' Fits a Mallows mixture model to ranking data.
#' 
#' Fits the Mallows mixture model to full ranking data, using
#' Kendall distance and an EM algorithm, for clustering permutations.
#' 
#' 
#' @param datas Matrix of fully-ranked data or permutations.
#' @param G Number of modes, 2 or greater.
#' @param weights Vector of frequencies of each permutation observed.
#' @param iter Maximum number of iterations.
#' @param tol Stopping precision.
#' @param logsumexp.trick Logical. Use logsumexp trick to compute log-likelihood.
#' @return See output of FormatOut.
#' @author Yunlong Jiao
#' @export
#' @references "Mixtures of distance-based models for ranking data". Thomas 
#' Brendan Murphy & Donal Martin. 1 April 2002. Computational Statistics & 
#' Data Analysis 41 (2003) 645-655.
#' @keywords cluster Mallows mixture
#' 

Mallows <-
  function(datas, G, weights = NULL, iter = 100, tol = 1e-3, logsumexp.trick = TRUE, iterin = iter, 
           exhkey = "_Exh", eqlamkey = "_Eqlam", 
           key = c("copelandMallows", "bruteMallows", "bordaMallows", "kernelMallows", "kernelMallows_Exh", 
                   "copelandMallows_Eqlam", "bruteMallows_Eqlam", "bordaMallows_Eqlam", "kernelMallows_Eqlam", "kernelMallows_Exh_Eqlam", 
                   "kernelGaussian", "kernelGaussian_Eqlam"))
  {
    if (is.null(weights)) {
      weights <- rep(1, nrow(datas))
    }
    
    key <- match.arg(key)
    exh <- grepl(exhkey, key)
    eqlam <- grepl(eqlamkey, key)
    if (grepl("kernelMallows", key)) {
      return(Mallows_kernel(datas = datas, G = G, iter = iter, weights = weights, eqlam = eqlam, tol = tol, iterin = iterin, exh = exh, logsumexp.trick = logsumexp.trick, key = key))
    } else if (grepl("kernelGaussian", key)) {
      return(gaussianmix(datas = datas, G = G, iter = iter, weights = weights, eqlam = eqlam, tol = tol, cov.constr = TRUE, key = key))
    }
    
    # Number of subjects.
    N <- nrow(datas)
    # Number of items being ranked.
    abils <- ncol(datas)
    
    if (grepl("bruteMallows", key)) {
      perm <- do.call("rbind", permn(abils))
    } else {
      perm <- NULL
    }
    
    dists.table <- DistanceDistribution(abils)
    # Initialize the p-value of membership in each cluster.
    p <- rep(1/G, G)
    # Initialize the modal sequences.
    R <- lapply(1:G, function(i) sample(abils))
    # Initialize the lambda values.
    lambda <- rep(runif(1, min = 0.1, max = 10), G)
    C.lam <-  unlist(lapply(lambda, 
                            function(i) C_lam(i, dists.table = dists.table, return.logC = logsumexp.trick)))
    
    # cat("Solving...\n")
    likelihood <- 0*(1:iter)
    infos <- KendallInfo(datas)
    all.dists.data <- AllKendall(datas, do.call("rbind", R), infos)
    i <- 1
    
    while (i <= iter) {
      # E Step
      z <- EStep(R, datas, p, lambda, G, N, C.lam, all.dists.data, use.logC = logsumexp.trick)
      
      # M Step
      R <- UpdateR(datas, weights * z, infos, perm, key = key)
      
      p <- UpdateP(z, weights)
      
      all.dists.data <- AllKendall(datas, do.call("rbind", R), infos)
      lambda <- UpdateLambda(r = datas, R = R, z = weights * z, G = G, all.dists.data = all.dists.data,
                             dists.table = dists.table, eqlam = eqlam)
      C.lam <- unlist(lapply(lambda, 
                             function(i) C_lam(i, dists.table = dists.table, return.logC = logsumexp.trick)))
      
      likelihood[i] <- Likelihood(z, p, C.lam, lambda, all.dists.data, weights, use.logC = logsumexp.trick)
      if (i > 2 && is.finite(likelihood[i]) && is.finite(likelihood[i-1]) && abs(likelihood[i] - likelihood[i-1]) < tol) {
        # message("Algorithm converged")
        likelihood[i:iter] <- likelihood[i]
        i <- iter
      }
      i <- i + 1
    }
    # cat("Formatting output")
    out <- FormatOut(R, p, lambda, z, datas, likelihood, all.dists.data, weights, key)
    return(out)
  }
