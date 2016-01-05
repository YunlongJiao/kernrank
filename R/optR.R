#' @importFrom stats nlm
#' 

optR <- function(rold, zz, lam, data.info, perm.info)
{
  negloglike_by_r <- function(r, zz, lam, data.info, perm.info){
    dg <- as.vector(distL2(r = perm.info, centers = matrix(r, nrow = 1)))
    ad <- as.vector(distL2(r = data.info, centers = matrix(r, nrow = 1)))
    ff <- -lam * dg
    logsum <- logsumexp(ff)
    negloglike <- lam * sum(zz * ad) + sum(zz) * logsum
    attr(negloglike, "gradient") <- 2 * lam * (sum(zz) * colSums(exp(ff-logsum)*perm.info) - colSums(zz*data.info))
    negloglike
  }
  
  res <- try(stats::nlm(f = negloglike_by_r, p = rold, zz = zz, lam = lam, data.info = data.info, perm.info = perm.info)$estimate, silent = TRUE)
  if (inherits(res, "try-error")) {
    warning("R not fully optimized due to stats::nlm error")
    return(rold)
  } else {
    return(res)
  }
}
