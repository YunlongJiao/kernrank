#' Kendall kernel for interleaving partial rankings
#' 
#' Calculates Kendall kernel between interleaving partial rankings in time O (klogk),
#' where ties (supposed few) are broken by adopting a convolution kernel
#' averaging compatible rankings without ties.
#' 
#' 
#' @param x Vector. 
#' If x is numeric, the rank vector converted from x indicate that larger values mean being preferred.
#' NAs replace unobserved values. 
#' @param y Same as x. 
#' @return Kendall kernel for interleaving partial rankings defined as kendall tau averaged over all compatible full ranking.
#' @author Yunlong Jiao
#' @export
#' @keywords Kendall Kernel PartialRanking
#' @references 
#' Jiao, Y., & Vert, J.-P. (2016). The Kendall and Mallows Kernels for Permutations. 2016. \href{https://hal.archives-ouvertes.fr/hal-01279273}{hal-01279273}
#' @examples 
#' x <- c(1.5, 0.1, NA, -4, NA)
#' y <- c(NA, NA, 0, 3, NA)
#' kendall_partial(x, y)
#' 

kendall_partial <- function(x, y)
{
  if (!is.vector(x) || !is.vector(y)) 
    stop("x and y must be vectors")
  if (any(is.na(x)) || any(is.na(y))) {
    FUN <- kendall_partial_inner
  } else {
    FUN <- kendall_corr
  }
  
  if (xd <- (anyDuplicated(x[is.finite(x)]) > 0)) {
    x <- dupgen(x)
  }
  if (yd <- (anyDuplicated(y[is.finite(y)]) > 0)) {
    y <- dupgen(y)
  }
  
  if (!xd && !yd) {
    return(FUN(x, y))
  } else if (xd && !yd) {
    return(mean(apply(x, 1, FUN, y)))
  } else if (!xd && yd) {
    return(mean(apply(y, 1, FUN, x)))
  } else {
    return(mean(apply(x, 1, function(v) apply(y, 1, FUN, v))))
  }
}
