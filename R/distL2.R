
distL2 <- function(r, centers, mc = 0.25)
{
  # @param r,centers Must be results of calling \code{\link{KendallInfo}}
  # @seealso \code{\link{KendallInfo}}
  # @return Squared euclidean distance in feature space multiplied by "\code{mc}"
  
  stopifnot(ncol(r)==ncol(centers))
  
  dists <- matrix(0, nrow = nrow(r), ncol = nrow(centers))
  for(i in 1:nrow(centers)){
    dists[ ,i] <- colSums((t(r) - centers[i, ])^2)
  }
  
  return(dists*mc)
}
