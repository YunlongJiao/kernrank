#' Common ranking aggregation methods for permutations
#' 
#' Used to update modal sequences of each cluster in the EM algorithm when fitting Mallows mixture models.
#' 
#' 
#' @param r A vector or a matrix of sequences in rows.
#' @param z A vector of weights/frequencies of observations or a matrix of probability of cluster membership 
#' for each sequence and each cluster. Set by default a constant vector of 1.
#' @param infos The result of KendallInfo(r). Optional for speeding up computation. 
#' @param perm A matrix of full permutations for brute-force search of optimal Kemeny consensus. Only effective for "key" equal to "brute".
#' @param key A character string indicating the ranking aggregation method to find centers.
#' \itemize{
#' \item \emph{borda} denotes the Borda count
#' \item \emph{copeland} denotes the Copeland's aggregated ranking
#' \item \emph{brute} denotes the optimal Kemeny consensus found by brute-force search
#' }
#' @return List of length 1 if "z" is a vector, or ncol(z) if "z" is a matrix, each entry being the modal sequence in each cluster.
#' @author Yunlong Jiao
#' @export 
#' @importFrom combinat permn
#' @keywords RankAggregation TotalRanking
#' @examples 
#' r <- do.call("rbind", list(1:5, 5:1, c(2,4,1,5,3)))
#' RankAggreg(r, key = "borda") # Borda count for sequences in "r"
#' 

RankAggreg <- function(r, z = NULL, infos = NULL, perm = NULL, key = c("borda", "copeland", "brute"))
{
  if (is.vector(r)) {
    r <- matrix(r, nrow = 1)
  } else if (!is.matrix(r)) {
    r <- as.matrix(r)
  }
  
  if (is.null(z)) {
    z <- rep(1, nrow(r))
  }
  if (is.vector(z)) {
    z <- matrix(z, ncol = 1)
  }
  
  if (grepl("brute", key)) {
    # DO NOT RUN for over-sized permutation !
    abils <- ncol(r)
    if (is.null(perm)) {
      perm <- do.call("rbind", combinat::permn(abils))
    }
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
    # result of the following heavily depends on the column order of combn results as ties are broken by FCFS-principle
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
    stop("Unable to pick center update function. Please specify ", sQuote("key"), " to be one of the implemented methods!")
  }
  return(R)
}
