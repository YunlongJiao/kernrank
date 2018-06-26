#' Weighted Kendall kernel for total rankings
#' 
#' Calculates weighted Kendall kernel between total rankings in time \code{O(nlogn)},
#' where ties (supposed few) are broken by adopting a convolution kernel
#' averaging compatible rankings without ties.
#' 
#' 
#' @param x,y Vector. 
#' If \code{x} is numeric, the rank vector converted from \code{x} indicate that larger values mean being preferred.
#' NAs are not allowed.
#' @param method Character. The method to perform weighted Kendall kernel. Choices include
#' \itemize{
#' \item \code{aken} denotes the Average Kendall kernel.
#' \item \code{ken} denotes the standard Kendall kernel.
#' \item \code{top} denotes the TOP-k Kendall kernel.
#' \item \code{add} denotes weighted Kendall kernel with ADDitive weights.
#' \item \code{mult} denotes weighted Kendall kernel with MULTiplicative weights.
#' }
#' @param k Integer. The parameter in top-\code{k} Kendall kernel.
#' Top-\code{k} implies ranks larger than k, where larger ranks mean being more preferred.
#' @param u Numeric vector. The parameter in additive or multiplicative weighted Kendall kernel.
#' @param normalized Logical. Whether to normalize the output kernel value.
#' The weighted Kendall kernel elaborated in Jiao and Vert (2018) corresponds to
#' the non-normalized version by setting \code{normalized=FALSE}.
#' @return Weighted Kendall kernel for total rankings,
#' where ties (supposed few) are broken by adopting a convolution kernel averaging compatible rankings without ties.
#' @author Yunlong Jiao
#' @export
#' @keywords Weighted Kendall Kernel TotalRanking
#' @references 
#' Yunlong Jiao, Jean-Philippe Vert. "The Weighted Kendall and High-order Kernels for Permutations." arXiv preprint arXiv:1802.08526, 2018. \href{https://arxiv.org/abs/1802.08526}{arXiv:1802.08526}
#' @examples 
#' x <- c(1.5, 0.1, 0, -4, 0)
#' y <- c(0, 0, 0, 3, 0)
#' 
#' # Average Kendall kernel
#' kendall_weight(x, y, method = "aken")
#' 
#' # Top-1 Kendall kernel is equiv to Standard Kendall kernel
#' kendall_weight(x, y, method = "ken")
#' kendall_weight(x, y, method = "top", k = 1)
#' 
#' # Additive/multiplicative weights with hyperbolic reduction factor
#' u <- 1/(length(x):1 + 1)
#' kendall_weight(x, y, method = "add", u = u)
#' kendall_weight(x, y, method = "mult", u = u)
#' 
#' # NOTE: In case of ties, standard Kendall kernel implemented by kendall_weight
#' # is NOT equiv to kendall_total due to different normalization technique!
#' # kendall_weight breaks ties by averaging over compatible full rankings:
#' kendall_weight(x, y, method = "ken", normalized = TRUE)
#' # kendall_total implements type-b of Kendall's tau that handles ties per se:
#' kendall_total(x, y)
#' 

kendall_weight <- function(x, y, 
                           method = c("aken", "ken", "top", "add", "mult"), 
                           k = NULL, u = NULL, normalized = TRUE)
{
  if (!is.vector(x) || !is.vector(y))
    stop("\"x\" and \"y\" must be vectors")
  
  if (any(is.na(x)) || any(is.na(y)))
    stop("\"x\" and \"y\" cannot contain any NAs!")
  
  # For ease of programming, convert R-readable obj to C-readable obj
  method <- match.arg(method)
  key <- switch(
    method,
    "ken" = 1,
    "aken" = 2,
    "top" = 3,
    "add" = 4,
    "mult" = 5
  )
  
  if (method == "top") {
    if (is.null(k))
      stop("\"k\" needs to be an integer for top-k Kendall kernel!")
    else
      k <- as.integer(k)
  } else {
    k <- 0L
  }
  
  if (method == "add" || method == "mult") {
    if (is.null(u) || length(u) != length(x))
      stop("\"u\" must be supplied as a numeric vector of the same length as \"x\"!")
    else
      u <- as.numeric(u)
  } else {
    u <- rep(0, length(x))
  }
  
  if (xd <- (anyDuplicated(x) > 0)) {
    x <- dupgen(x) # duplicates removed
  }
  if (yd <- (anyDuplicated(y) > 0)) {
    y <- dupgen(y) # duplicates removed
  }
  
  if (!xd && !yd) {
    return(kendall_weight_inner(x, y, key, k, u, normalized))
  } else if (xd && !yd) {
    return(mean(apply(x, 1, kendall_weight_inner, y, key, k, u, normalized)))
  } else if (!xd && yd) {
    return(mean(apply(y, 1, kendall_weight_inner, x, key, k, u, normalized)))
  } else {
    return(mean(apply(x, 1, function(v) apply(y, 1, kendall_weight_inner, v, key, k, u, normalized))))
  }
}
