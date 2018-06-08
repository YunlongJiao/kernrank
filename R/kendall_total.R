#' Kendall kernel for total rankings
#' 
#' Calculates Kendall kernel between total rankings in time \code{O(nlogn)},
#' where ties are dealt with type-b (soft version) of kendall kernel
#' 
#' 
#' @param x A vector or a matrix of m1 sequences in rows and orders of n items in cols.
#' If \code{x} is numeric, the rank vector converted from \code{x} indicate that larger values mean being preferred.
#' NAs are not allowed.
#' @param y Same as \code{x}. By default "\code{y}" is set equal to "\code{x}".
#' @return Kendall kernel for total rankings defined as type-b (soft version) of kendall kernel.
#' @export
#' @keywords Kendall Kernel TotalRanking
#' @references 
#' Kendall rank correlation coefficient: \url{https://en.wikipedia.org/wiki/Kendall_rank_correlation_coefficient}
#' @references 
#' Yunlong Jiao, Jean-Philippe Vert. "The Kendall and Mallows Kernels for Permutations." IEEE Transactions on Pattern Analysis and Machine Intelligence (TPAMI), vol. 40, no. 7, pp. 1755-1769, 2018. \href{https://doi.org/10.1109/TPAMI.2017.2719680}{DOI:10.1109/TPAMI.2017.2719680}
#' @examples 
#' x <- c(1.5, 0.1, 0, -4, 0)
#' y <- c(0, 0, 0, 3, 0)
#' kendall_total(x, y)
#' 

kendall_total <- function(x, y = NULL)
{
  if (is.null(y) && is.vector(x)) {
    kmat <- kendall_corr(x, x) # kmat is number
  } else if (is.null(y) && !is.vector(x)) {
    x <- as.matrix(x)
    kmat <- kendall_corr(t(x)) # kmat is matrix
    dimnames(kmat) <- list(rownames(x), rownames(x))
  } else if (is.vector(y) && is.vector(x)) {
    kmat <- kendall_corr(x, y) # kmat is number
  } else {
    # Convert both x and y to matrix
    if (is.vector(x)) {
      x <- matrix(x, nrow = 1)
    } else if (!is.matrix(x)) {
      x <- as.matrix(x)
    }
    if (is.vector(y)) {
      y <- matrix(y, nrow = 1)
    } else if (!is.matrix(y)) {
      y <- as.matrix(y)
    }
    # Generate kernel matrix
    if (ncol(x) != ncol(y)) 
      stop("\"x\" and \"y\" need to have the same number of ranked items!")
    kmat <- kernmat(kf = kendall_corr, mat1 = x, mat2 = y) # kmat is matrix
  }
  return(kmat)
}
