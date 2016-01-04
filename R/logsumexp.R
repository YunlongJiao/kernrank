#' Computing Log-Sum-Exp with a common trick
#' 
#' Computes log-sum-exp of a series of (typically large in absolute value) numbers with a trick
#' 
#' 
#' @param x A matrix or vector of numerics.
#' @param byrow Logical. Computes by rows if a matrix is provided.
#' @param bycol Logical. Computes by cols if a matrix is provided.
#' @return Log-Sum-Exp of the numbers in vector or row-/column-wise Log-Sum-Exp of the numbers in matrix.
#' @author Yunlong Jiao
#' @export
#' @references https://hips.seas.harvard.edu/blog/2013/01/09/computing-log-sum-exp/
#' @examples
#' 
#' ## Not run:
#' 
#' x <- c(-1000, -999, -1000)
#' logsumexp(x)
#' 
#' ## End(Not run)

logsumexp <- function(x, byrow = TRUE, bycol = !byrow)
{
  if (is.matrix(x)) {
    if (byrow && !bycol) {
      a1 <- apply(x, 1, max)
      res <- a1 + log(rowSums(exp(x - a1)))
    } else if (bycol) {
      a1 <- apply(x, 2, max)
      res <- a1 + log(rowSums(exp(t(x) - a1)))
    }
  }
  if (is.vector(x) || (is.matrix(x) && !byrow && !bycol)) {
    a1 <- max(x)
    res <- a1 + log(sum(exp(x - a1)))
  }
  res
}
