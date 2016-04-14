#' All Kendall's distances between two sets of total rankings or real-valued vectors
#' 
#' Calculates all of the Kendall's distances between two sets of total 
#' rankings or real-valued vectors. 
#' 
#' 
#' @param r A vector or a matrix of m1 sequences in rows and orders of n items in cols.
#' @param seqs Another vector or a matrix of m2 sequences in rows and orders of n items in cols. By default "seqs" is set equal to "r".
#' @param data.info Optional argument giving the Kendall embedding of "r", that is the result of KendallInfo(r), 
#' to facilitate computing Kendall's difference for "r" to "seqs" without exploring the kernel trick.
#' @param use.kernel.trick Logical indicating whether the kernel trick is explored. 
#' This is particularly interesting when the number of items to be ranked is high 
#' and pcaPP::cor.fk() is available. 
#' By default (set FALSE), Kendall embedding is explicitly computed; 
#' otherwise kernel trick is explored. 
#' @param kmat Kendall kernel matrix of dimension m1 x m2, correlation type correponding to "type". If given, kernel trick is explored directly.
#' @param type A character string indicating the type of Kendall correlation for "kmat".
#' @param mc A normalization constant default to 0.25 such that output normalized squared Euclidean distance in the feature space induced by Kendall embedding amounts exactly to Kendall distances.
#' @return A matrix of dimension m1 x m2 where entry [i,j] is the distance from sequence "i"
#' in "r" to sequence "j" in "seqs".
#' @note Kernel trick is explored in the sense that "r" and "seq" are only used for checking dimensions and getting attributes but not used explicitly to compute the distance.
#' This is particularly interesting when data is high-dimensional in constrast to rather few observations (m1,m2>>n). 
#' Option "use.kernel.trick" set TRUE or FALSE may give slightly different results due to computation precision.
#' @author Yunlong Jiao
#' @keywords Kendall Distance TotalRanking
#' @export
#' @references 
#' Kendall rank correlation coefficient: \url{https://en.wikipedia.org/wiki/Kendall_rank_correlation_coefficient}
#' @references 
#' Jiao, Y., & Vert, J.-P. (2016). The Kendall and Mallows Kernels for Permutations. 2016. \href{https://hal.archives-ouvertes.fr/hal-01279273}{hal-01279273}
#' @examples
#' 
#' #### Ex 1: compute Kendall distance matrix and Mallows kernel matrix
#' data1 <- do.call("rbind", list(1:5, 5:1, c(3, 2, 1, 4, 5)))
#' data2 <- do.call("rbind", list(1:5, 5:1))
#' 
#' # Kendall distance matrix
#' s.K.d.mat <- AllKendall(data1, data2)
#' 
#' # Mallows kernel matrix with dispersion parameter lambda
#' lambda <- 0.1
#' M.k.kmat <- exp(-lambda * s.K.d.mat)
#' 
#' #### Ex 2: why kernel trick?
#' r <- lapply(1:20, function(i) sample.int(1000, replace = TRUE))
#' r <- do.call('rbind', r)
#' dim(r)
#' 
#' # I) without kernel trick
#' pt <- proc.time()
#' dmat1 <- AllKendall(r, use.kernel.trick = FALSE)
#' proc.time() - pt
#' 
#' # II) with kernel trick (should be much faster in this setting)
#' require(pcaPP)
#' pt <- proc.time()
#' dmat2 <- AllKendall(r, use.kernel.trick = TRUE)
#' proc.time() - pt
#' 
#' # NOTE: dmat1 and dmat2 may return slightly different values due to computation precision
#' isTRUE(all.equal(dmat1, dmat2, check.attributes = FALSE)) # may sometimes output FALSE
#' isTRUE(max(abs(dmat1 - dmat2)) < 1e-6) # always output TRUE
#' 

AllKendall <- function(r, seqs = NULL, data.info = NULL, use.kernel.trick = FALSE, 
                       kmat = NULL, type = c("type-b", "type-a"), mc = 0.25)
{
  type <- match.arg(type)
  
  if (is.vector(r)) {
    r <- matrix(r, nrow = 1)
  } else if (!is.matrix(r)) {
    r <- as.matrix(r)
  }
  
  if (!use.kernel.trick && is.null(kmat)) {
    # CASE I compute explicitly all Kendall embeddings
    if (is.null(data.info)) {
      data.info <- KendallInfo(r)
    }
    if (is.null(seqs)) {
      dists <- distL2(r = data.info, centers = data.info, mc = mc)
    } else {
      stopifnot(ncol(r) == ncol(seqs))
      seqs.info <- KendallInfo(seqs)
      dists <- distL2(r = data.info, centers = seqs.info, mc = mc)
      return(dists)
    }
  } else {
    if (is.null(kmat)) {
      # CASE I-1 explore kernel trick in high-dimensional setting
      type <- 'type-b'
      kmat <- kendall_total(x = r, y = seqs)
    }
    # CASE II if kmat is given then kernel trick is explored directly no matter what
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
}
