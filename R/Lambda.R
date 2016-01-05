
Lambda <- function(lambda, rhs, dists = NULL, dists.table = NULL)
{
  if (!is.null(dists)) {
    ps <- -lambda * dists
    logsum <- -logsumexp(ps, byrow = FALSE, bycol = FALSE)
    lhs <- sum(dists * exp(logsum - lambda * dists))
  } else if (!is.null(dists.table)) {
    lhs <- C_lam(lambda, 
                 dists.table = dists.table, return.logC = FALSE) * sum(as.numeric(names(dists.table))*dists.table*exp(-lambda*as.numeric(names(dists.table))))
  } else {
    stop("Unable to calculate Lambda function due to both dists and dists.table being NULL")
  }
  lhs - rhs
}
