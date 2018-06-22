#' @importFrom stats cor
#' 

kendall_corr <- function(x, y)
{
  # Kendall correlation function for two vectors
  if (any(is.na(x)) || any(is.na(y)))
    stop("\"x\" and \"y\" cannot contain any NAs!")
  
  if (requireNamespace("pcaPP", quietly = TRUE)) {
    pcaPP::cor.fk(x = x, y = y)
  } else {
    stats::cor(x = x, y = y, use = "everything", method = "kendall")
  }
}
