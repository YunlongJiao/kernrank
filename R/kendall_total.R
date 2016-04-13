#' Kendall kernel for total rankings
#' 
#' Calculates Kendall kernel between total rankings in time O (nlogn),
#' where ties are dealt with type-b (soft version) of kendall kernel
#' 
#' 
#' @param x Numeric vector. A rank vector is converted from x then,
#' where larger values indicate being preferred and NA not allowed.
#' @param y Same as x.
#' @return Kendall kernel for total rankings defined as type-a (soft version) of kendall kernel.
#' @importFrom pcaPP cor.fk
#' @export
#' @note The implementation is directly taken from pcaPP::cor.fk()
#' @keywords kendall kernel total ranking
#' @references 
#' Kendall rank correlation coefficient: \url{https://en.wikipedia.org/wiki/Kendall_rank_correlation_coefficient}
#' @references 
#' Jiao, Y., & Vert, J. P. (2015). The Kendall and Mallows kernels for permutations. In Proceedings of the 32nd International Conference on Machine Learning (ICML-15) (pp. 1935-1944).
#' @examples 
#' x <- c(1.5, 0.1, 0, -4, 0)
#' y <- c(0, 0, 0, 3, 0)
#' kendall_total(x, y)
#' 

kendall_total <- function(x, y)
{
  if (any(is.na(x)) || any(is.na(y))) {
    stop("NA not allowed in kendall kernel for total rankings")
  } else {
    FUN <- pcaPP::cor.fk
  }
  
  return(FUN(x, y))
}
