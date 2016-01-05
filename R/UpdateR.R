#' Common ranking aggregation methods.
#' 
#' Used in EM algorithm for fitting Mallows model to update modal sequences 
#' in each cluster.
#' 
#' @param r Matrix of sequences being clustered. Permutations in rows.
#' @param z Vector of weights or matrix of probability of cluster membership 
#' for each sequence and each cluster.
#' @param infos The KendallInfo matrix for "r". Optional for speeding up computation. 
#' @param perm Full list of permutations for brute-force search of optimal center.
#' @param key Character. Ranking aggregation method to define centers.
#' @return List. Consensus ranking or new cluster centers for each cluster.
#' @author Yunlong Jiao
#' @export
#' @keywords cluster center
#' @examples 
#' ## Not run:
#' 
#' r <- do.call("rbind", list(1:5, 5:1, c(2,4,1,5,3)))
#' UpdateR(r, key = "borda")
#' 
#' ## End(Not run)

UpdateR <- function(r, z = NULL, infos = NULL, perm = NULL, key = c("borda", "copeland", "brute"))
{
  if (is.null(z)) {
    z <- rep(1, nrow(r))
  }
  
  if (is.vector(z)) {
    z <- matrix(z, ncol = 1)
  }
  
  if (grepl("brute", key)) {
    # DO NOT RUN for over-sized permutation / rewrite !
    dists <- AllKendall(r = r, seqs = perm, data.info = infos)
    R <- lapply(1:ncol(z), function(j){
      idx <- which.min(colSums(dists * z[ ,j]))
      perm[idx, ]
    })
  } else if (grepl("borda", key)) {
    R <- lapply(1:ncol(z), function(j){
      cent <- apply(r * z[ ,j], 2, mean)
      rank(cent, ties.method = "first")
    })
  } else if (grepl("copeland", key)) {
    if (is.null(infos)) {
      infos <- KendallInfo(r)
    }
    abils <- ncol(r)
    # the following strictly dependent on the column order of combn results
    num <- c(0, cumsum((abils-1):1))
    R <- lapply(1:ncol(z), function(j){
      cent.info <- sign(colSums(infos * z[ ,j]))
      cent <- 0*(1:abils)
      for (i in 1:(abils-1)) {
        cent[i] <- cent[i] + sum(cent.info[(num[i]+1):(num[i+1])])
      }
      for (i in 2:abils) {
        cent[i] <- cent[i] - sum(cent.info[num[2:i]-(abils-i)])
      }
      rank(cent, ties.method = "first")
    })
  } else {
    stop("Unable to pick UpdateR function. Please specify ", sQuote("key"), " to be one of the following : borda, copeland, brute!")
  }
  return(R)
}
