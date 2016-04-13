#' The E-step of the EM algorithm
#' 
#' Assigns each ranking the probability that it belongs to each cluster, given
#' current parameters.
#' 
#' 
#' @param R List of length "G". Current cluster modal sequences.
#' @param r A matrix with "N" sequences of full rankings in rows.
#' @param p Numeric vector. The proportion probability of each cluster.
#' @param lambda Numeric vector. The lambda parameters from Mallows model for each cluster.
#' @param G Number of clusters, equal to length(R).
#' @param N Number of rows in the data.
#' @param C Vector of normalizing coefficients for the clusters.
#' @param all.dists For efficiency, provide all of the Kendall distances
#' between each sequence and each cluster mode, that is the result of AllKendall(r,do.call('rbind',R)).
#' @param use.logC Logical. Whether or not to use log-sum-exp trick.
#' @return Matrix of dimension N x G where output[i,j] represents the current fuzzy assignment membership probability, 
#' namely the updated probability that subject "i" belongs to cluster "j".
#' @author Yunlong Jiao
#' @references 
#' Murphy, T. B., & Martin, D. (2003). Mixtures of distance-based models for ranking data. Computational statistics & data analysis, 41(3), 645-655.
#' @keywords expectation maximization
#' 

EStep <-
function(R, r, p, lambda, G, N, C, all.dists, use.logC = FALSE) {
  if (use.logC) {
    fs <- exp(C-lambda*t(all.dists))
  } else {
    fs <- C*exp(-lambda*t(all.dists))
  }
  z <- t(p*fs)
  
  z <- z/rowSums(z)
  return(z)
}
