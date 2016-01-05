#' Kendall kernel for top-k rankings
#' 
#' Calculates Kendall kernel between top-k rankings in time O (klogk),
#' where ties (supposed few) are broken by adopting a convolution kernel
#' averaging compatible rankings without ties.
#' 
#' 
#' @param x Numeric vector. A rank vector is converted from x then,
#' where larger values indicate being preferred and NA replace unobserved values. 
#' @param y Same as x.
#' @return Kendall kernel for top-k rankings defined as kendall tau averaged over all compatible full ranking.
#' @importFrom pcaPP cor.fk
#' @author Yunlong Jiao
#' @export
#' @keywords kendall kernel top-k ranking
#' @examples 
#' ## Not run:
#' 
#' x <- c(1.5, 0.1, NA, -4, NA)
#' y <- c(NA, NA, 0, 3, NA)
#' kendall_top(x, y)
#' 
#' ## End(Not run)

kendall_top <- function(x, y)
{
  
  if (any(is.na(x)) || any(is.na(y))) {
    FUN <- kendall_top_inner
  } else {
    FUN <- pcaPP::cor.fk
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