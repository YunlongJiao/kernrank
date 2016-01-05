#' The E-step of the EM algorithm.
#' 
#' Assigns each ranking the probability that it belongs to each cluster, given
#' current parameters.
#' 
#' 
#' @param R List. Current cluster modal sequences.
#' @param r The data of full rankings.
#' @param p The proportion of the data currently assigned to each cluster.
#' @param lambda The lambda parameters from Mallows model for each cluster.
#' @param G Number of clusters, length(R).
#' @param N Number of rows in the data.
#' @param C Vector of normalizing coefficients for the clusters.
#' @param all.dists For efficiency, provide all of the Kendall distances
#' between each sequence and each cluster mode.
#' @return Matrix where output[i, j] represents the current probability that
#' subject "i" belongs to cluster "j".
#' @author Yunlong Jiao
#' @references 
#' "Mixtures of distance-based models for ranking data". Thomas 
#' Brendan Murphy & Donal Martin. 1 April 2002. Computational Statistics & 
#' Data Analysis 41 (2003) 645-655.
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
