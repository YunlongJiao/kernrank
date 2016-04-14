#' Kendall kernel for total rankings
#' 
#' Calculates Kendall kernel between total rankings in time O (nlogn),
#' where ties are dealt with type-b (soft version) of kendall kernel
#' 
#' 
#' @param x A vector or a matrix of m1 sequences in rows and orders of n items in cols.
#' If x is numeric, the rank vector converted from x indicate that larger values mean being preferred.
#' NAs are not allowed.
#' @param y Same as x. By default "y" is set equal to "x".
#' @return Kendall kernel for total rankings defined as type-b (soft version) of kendall kernel.
#' @export
#' @keywords Kendall Kernel TotalRanking
#' @references 
#' Kendall rank correlation coefficient: \url{https://en.wikipedia.org/wiki/Kendall_rank_correlation_coefficient}
#' @references 
#' Jiao, Y., & Vert, J.-P. (2016). The Kendall and Mallows Kernels for Permutations. 2016. \href{https://hal.archives-ouvertes.fr/hal-01279273}{hal-01279273}
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
    # convert both x and y to matrix
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
    # generate kernel matrix
    if (ncol(x) != ncol(y)) 
      stop("x and y need to have the same number of ranked items")
    kmat <- kernmat(kf = kendall_corr, mat1 = x, mat2 = y) # kmat is matrix
  }
  return(kmat)
}
