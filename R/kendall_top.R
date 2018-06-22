#' Kendall kernel for top-k rankings
#' 
#' Calculates Kendall kernel between top-k rankings in time \code{O(klogk)},
#' where ties (supposed few) are broken by adopting a convolution kernel
#' averaging compatible rankings without ties.
#' 
#' 
#' @param x,y Vector. 
#' If \code{x} is numeric, the rank vector converted from \code{x} indicate that larger values mean being preferred.
#' NAs replace unobserved values.
#' @return Kendall kernel for top-k rankings,
#' where ties (supposed few) are broken by adopting a convolution kernel averaging compatible rankings without ties.
#' @author Yunlong Jiao
#' @export
#' @keywords Kendall Kernel TopkRanking
#' @references 
#' Yunlong Jiao, Jean-Philippe Vert. "The Kendall and Mallows Kernels for Permutations." IEEE Transactions on Pattern Analysis and Machine Intelligence (TPAMI), vol. 40, no. 7, pp. 1755-1769, 2018. \href{https://doi.org/10.1109/TPAMI.2017.2719680}{DOI:10.1109/TPAMI.2017.2719680}
#' @examples 
#' x <- c(1.5, 0.1, NA, -4, NA)
#' y <- c(NA, NA, 0, 3, NA)
#' kendall_top(x, y)
#' 

kendall_top <- function(x, y)
{
  if (!is.vector(x) || !is.vector(y)) 
    stop("\"x\" and \"y\" must be vectors")
  if (any(is.na(x)) || any(is.na(y))) {
    FUN <- kendall_top_inner
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
