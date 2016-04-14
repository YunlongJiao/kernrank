
EStep <-
function(R, r, p, lambda, G, N, C, all.dists, use.logC = FALSE) {
  if (use.logC) {
    fs <- exp(C-lambda*t(all.dists))
  } else {
    fs <- C*exp(-lambda*t(all.dists))
  }
  z <- t(p*fs)
  
  z <- z/rowSums(z)
  return(z)
}
