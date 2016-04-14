
Likelihood <- function(z, p, C.lam, lambda, all.dists.data, weights, use.logC = FALSE)
{
  if (use.logC) {
    like <- log(p) + C.lam - lambda * t(all.dists.data)
    like <- LogSumExp(like, bycol = TRUE)
  } else {
    like <- log(colSums(C.lam * p * exp(-lambda * t(all.dists.data))))
  }
  sum(weights * like)
}
