#' Update modal sequences in each cluster.
#' 
#' Maximizes the likelihood of the data by updating the cluster centers of the
#' model.
#' 
#' @param r Matrix of sequences being clustered.
#' @param z Probability of cluster membership for each sequence and each
#' cluster.
#' @param infos The KendallInfo matrix for "r".
#' @param perm Full list of permutations for brute search.
#' @param key Character. Rank aggregation method to define centers.
#' @return New cluster centers for each cluster.
#' @author Yunlong Jiao
#' @references "Mixtures of distance-based models for ranking data". Thomas 
#' Brendan Murphy & Donal Martin. 1 April 2002. Computational Statistics & 
#' Data Analysis 41 (2003) 645-655.
#' @keywords cluster center
#' 
#' 

UpdateR <- function(r, z, infos = NULL, perm = NULL, key)
{
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
    stop("Unable to pick UpdateR function")
  }
  return(R)
}
