
UpdateP <- function(z, weights)
{
  colSums(weights * z)/sum(weights)
}
