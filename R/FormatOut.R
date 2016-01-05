
FormatOut <- function(R, p, lambda, z, datas, likelihood, all.dists.data, weights, key)
{
  # Find which cluster the individual belongs to.
  clust <- apply(z, 1, which.max)
  
  datas <- data.frame(datas, clust = clust, weights = weights, pvals = z, dists = all.dists.data)
  out <- list(key = key, R = R, p = p, lambda = lambda, 
              datas = datas, min.like = likelihood)
  return(out)
}
