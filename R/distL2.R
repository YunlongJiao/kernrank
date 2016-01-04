#' Calculates squared Euclidean distances between two different 
#' sets of vectors.
#' 
#' Amounts to all Kendall's distances when two vectors equal
#' two sets of score-based pairwise relative ordering, 
#' i.e., KendallInfo.
#' 

#' 
#' @param r One set of sequences.
#' @param centers Another set of sequences.
#' @param mc A constant multiplied to the distance. Default to 0.25 such that
#' output amounts to all Kendall's distances when two vectors equal
#' two sets of pairwise relative ordering.
#' @return Matrix where output[i, j] represents the squred euclidean distance 
#' from sequence "i" in "r" to sequence "j" in "seqs", multiplied by mc.
#' @author Yunlong Jiao
#' @keywords Squared Euclidean distance

distL2 <- function(r, centers, mc = 0.25)
{
  stopifnot(ncol(r)==ncol(centers))
  
  dists <- matrix(0, nrow = nrow(r), ncol = nrow(centers))
  for(i in 1:nrow(centers)){
    dists[ ,i] <- colSums((t(r) - centers[i, ])^2)
  }
  
  return(dists*mc)
}
