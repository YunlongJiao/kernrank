#' @importFrom stats cor
#' 

kendall_corr <- function(x, y = NULL)
{
  # Kendall correlation function
  
  if (requireNamespace("pcaPP", quietly = TRUE)) {
    pcaPP::cor.fk(x = x, y = y)
  } else {
    stats::cor(x = x, y = y, use = "everything", method = "kendall")
  }
}
