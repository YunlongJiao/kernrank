#' Computing Log-Sum-Exp with a common trick
#' 
#' Computes log-sum-exp of a series of (typically large in absolute value) numbers 
#' with a more accurate computational trick typically useful for small values.
#' 
#' 
#' @param x A vector or a matrix of numerics (typicall very small).
#' @param byrow Logical. Computes by rows if a matrix "\code{x}" is provided.
#' @param bycol Logical. Computes by cols if a matrix "\code{x}" is provided.
#' @return Log-Sum-Exp of the numbers in the vector "\code{x}", that is log(sum(exp(x))), 
#' or row-/column-wise Log-Sum-Exp of the numbers in a matrix "\code{x}".
#' @author Yunlong Jiao
#' @export
#' @references 
#' Computing Log-Sum-Exp: \url{https://hips.seas.harvard.edu/blog/2013/01/09/computing-log-sum-exp/}
#' @examples
#' x <- c(-1000, -999, -1000)
#' LogSumExp(x)
#' 

LogSumExp <- function(x, byrow = TRUE, bycol = !byrow)
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
