
C_lam <- function(lambda, dists.table = NULL, dists.to.Rg = NULL, return.logC = FALSE)
{
  if (!is.null(dists.table)) {
    # C_lam only works with single lambda value if dists.table given
    stopifnot(length(lambda) == 1)
    
    ps <- dists.table*(exp(-lambda*as.numeric(names(dists.table))))
    if (return.logC) {
      logC <- -log(sum(ps))
      return(logC)
    } else {
      C <- 1/sum(ps)
      return(C)
    }
  } else if (!is.null(dists.to.Rg)) {
    # C_lam only works with vectorized lambda if dists.to.Rg given
    stopifnot(length(lambda) > 1)
    
    ps <- -lambda * t(dists.to.Rg)
    if (return.logC) {
      logC <- -LogSumExp(ps, byrow = TRUE)
      return(logC)
    } else {
      C <- 1/rowSums(exp(ps))
      return(C)
    }
  } else {
    stop("unable to calculate normalization constant C_lam, at least either given dists.table or dists.to.Rg")
  }
}
