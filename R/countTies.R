
countTies <- function(x)
{
  # @param x Must be a matrix!
  # @return A vector representing number of tied pairs in each rows of "x".
  
  if (is.vector(x)) {
    x <- matrix(x, nrow = 1)
  } else if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  
  n1 <- apply(x, 1, function(u){
    tab <- table(u)
    sum(tab * (tab-1) / 2)
  })
  
  return(n1)
}
