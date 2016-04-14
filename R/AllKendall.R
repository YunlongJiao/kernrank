#' All Kendall's distances between two sets of total rankings or real-valued vectors
#' 
#' Calculates all of the Kendall's distances between two sets of total 
#' rankings or real-valued vectors.
#' 
#' @param r A vector or a matrix of m1 sequences in rows.
#' @param seqs Another vector or a matrix of m2 sequences in rows, by default set to equal to "r".
#' @param data.info Optional argument giving the Kendall embedding of "r", that is the result of KendallInfo(r), 
#' to facilitate computing Kendall's difference for "r" to "seqs".
#' @param high.dim Logical indicating whether the number of items to be ranked 
#' is high. By default (set FALSE), Kendall embedding is explicitly computed; 
#' otherwise kernel trick is explored.
#' @return A matrix of dimension m1 x m2 where entry [i,j] is the distance from sequence "i"
#' in "r" to sequence "j" in "seqs".
#' @note Option "high.dim" set TRUE or FALSE may give slightly different results due to computation precision.
#' @author Yunlong Jiao
#' @keywords Kendall distance
#' @importFrom pcaPP cor.fk
#' @importFrom kernlab kernelMatrix
#' @export
#' @references 
#' Kendall rank correlation coefficient: \url{https://en.wikipedia.org/wiki/Kendall_rank_correlation_coefficient}
#' @references 
#' Jiao, Y., & Vert, J.-P. (2016). The Kendall and Mallows Kernels for Permutations. 2016. \href{https://hal.archives-ouvertes.fr/hal-01279273}{hal-01279273}
#' @examples
#' # Ex 1: compute Kendall distance matrix and Mallows kernel matrix
#' 
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
#' # Ex 2: why kernel trick?
#' 
#' r <- lapply(1:20, function(i) sample.int(1000))
#' r <- do.call('rbind', r)
#' dim(r)
#' 
#' # I) without kernel trick
#' pt <- proc.time()
#' dmat1 <- AllKendall(r, high.dim = FALSE)
#' proc.time() - pt
#' 
#' # II) with kernel trick
#' pt <- proc.time()
#' dmat2 <- AllKendall(r, high.dim = TRUE)
#' proc.time() - pt
#' 
#' # NOTE: II should be around must faster than I in this setting
#' # NOTE: dmat1 and dmat2 may return slightly different values due to computation precision
#' 

AllKendall <- function(r, seqs = NULL, data.info = NULL, high.dim = FALSE)
{
  if (is.vector(r)) {
    r <- matrix(r, nrow = 1)
  } else if (!is.matrix(r)) {
    r <- as.matrix(r)
  }
  
  if (high.dim) {
    if (is.null(seqs)) {
      kmat <- pcaPP::cor.fk(t(r)) # type-b (soft version) of kendall corr with const diag equal to 1
    } else {
      kf <- pcaPP::cor.fk
      class(kf) <- 'kernel'
      kmat <- kernlab::kernelMatrix(kernel = kf, x = r, y = seqs)
    }
    dists <- distL2mat(kmat = kmat, r = r, seqs = seqs, type = "type-b")
  } else {
    if (is.null(data.info)) {
      data.info <- KendallInfo(r)
    }
    if (is.null(seqs)) {
      dists <- distL2(r = data.info, centers = data.info)
    } else {
      seqs.info <- KendallInfo(seqs)
      dists <- distL2(r = data.info, centers = seqs.info)
    }
  }
  
  return(dists)
}
