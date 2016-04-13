#' Convert a Kendall kernel matrix to a Kendall distance matrix
#' 
#' Normalized squared Euclidean distances in feature space recovers 
#' exactly the Kendall distance. 
#' 
#' @param kmat Kendall kernel matrix of dimension m1 x m2, correlation type correponding to "type".
#' @param r,seq Matrices of m1,m2 sequences in rows with n ranked items in cols, by default "seq" set equal to "r".
#' @param type A character string indicating the type of Kendall correlation.
#' @param mc A normalization constant default to 0.25 such that output squared Euclidean distance amounts exactly to Kendall distances.
#' @return A matrix of dimension m1 x m2 where entry [i,j] is the Kendall distance from sequence "i"
#' in "r" to sequence "j" in "seqs".
#' @note Kernel trick is explored in the sense that "r" and "seq" are only used for checking dimensions and getting attributes but not used explicitly to compute the distance.
#' This is particularly interesting when data is high-dimensional in constrast to rather few observations (m1,m2>>n).
#' @references 
#' Kendall rank correlation coefficient: \url{https://en.wikipedia.org/wiki/Kendall_rank_correlation_coefficient}
#' @references 
#' Jiao, Y., & Vert, J. P. (2015). The Kendall and Mallows kernels for permutations. In Proceedings of the 32nd International Conference on Machine Learning (ICML-15) (pp. 1935-1944).
#' @examples 
#' r <- lapply(1:20, function(i) sample.int(1000))
#' r <- do.call('rbind', r)
#' 
#' # Kendall kernel matrix first
#' kmat <- pcaPP::cor.fk(t(r))
#' # then Kendall distance matrix
#' dmat <- distL2mat(kmat, r, type = "type-b")
#' 
#' # notably the above is equivalent to directly calling AllKendall() with kernel trick
#' dmat <- AllKendall(r, high.dim = TRUE)
#' 

distL2mat <- function(kmat, r, seqs = NULL, type = c("type-b", "type-a"), mc = 0.25)
{
  type <- match.arg(type)
  
  if (is.null(seqs)) {
    seqs <- r
  }
  stopifnot(ncol(r) == ncol(seqs))
  n0 <- choose(ncol(r), 2)
  if (type == "type-a") {
    v1 <- rep(sqrt(n0), nrow(r))
    v2 <- rep(sqrt(n0), nrow(seqs))
  } else if (type == "type-b") {
    v1 <- sqrt(n0 - countTies(r))
    v2 <- sqrt(n0 - countTies(seqs))
  } else {
    stop("undefined type for kernel matrix!")
  }
  
  stopifnot(nrow(kmat) == length(v1))
  stopifnot(ncol(kmat) == length(v2))
  
  kmat <- sweep(
    sweep(kmat, 1, v1, "*"), 
    2, v2, "*")
  dmat <- sweep(
    sweep(-2*kmat, 1, v1*v1, "+"), 
    2, v2*v2, "+")
  return(mc*dmat)
}
