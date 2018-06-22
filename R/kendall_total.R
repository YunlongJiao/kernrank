#' Kendall kernel for total rankings
#' 
#' Calculates Kendall kernel between total rankings in time \code{O(nlogn)}, 
#' where ties are dealt with by type-b soft version of Kendall tau correlation.
#' 
#' 
#' @param x,y Vector.
#' If \code{x} is numeric, the rank vector converted from \code{x} indicate that larger values mean being preferred.
#' NAs are not allowed.
#' @return Kendall kernel for total rankings,
#' where ties are dealt with by type-b soft version of Kendall tau correlation.
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

kendall_total <- function(x, y)
{
  if (!is.vector(x) || !is.vector(y)) 
    stop("\"x\" and \"y\" must be vectors")
  
  if (any(is.na(x)) || any(is.na(y)))
    stop("\"x\" and \"y\" cannot contain any NAs!")
  
  return(kendall_corr(x, y))
}
