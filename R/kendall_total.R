#' Kendall kernel for total rankings
#' 
#' Calculates Kendall kernel between total rankings in time O (nlogn),
#' where ties are dealt with type-a (soft version) of kendall kernel
#' 
#' 
#' @param x Numeric vector. A rank vector is converted from x then,
#' where larger values indicate being preferred and NA not allowed.
#' @param y Same as x.
#' @return Kendall kernel for total rankings defined as type-a (soft version) of kendall kernel.
#' @importFrom pcaPP cor.fk
#' @export
#' @note The implementation is directly taken as pcaPP::cor.fk
#' @keywords kendall kernel total ranking
#' @examples 
#' ## Not run:
#' 
#' x <- c(1.5, 0.1, 0, -4, 0)
#' y <- c(0, 0, 0, 3, 0)
#' kendall_total(x, y)
#' 
#' ## End(Not run)

kendall_total <- function(x, y)
{
  if (any(is.na(x)) || any(is.na(y))) {
    stop("NA not allowed in kendall kernel for total rankings")
  } else {
    FUN <- pcaPP::cor.fk
  }
  
  return(FUN(x, y))
}