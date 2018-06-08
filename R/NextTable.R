
NextTable <- function(last.table, N.last)
{
  # @note Taken directly from R CRAN package \href{https://cran.r-project.org/web/packages/RMallow/index.html}{RMallow}.
  
  len <- (N.last + 1)*(N.last)/2 + 1
  # Length of vector, minus the stuff that is easy to fill in.
  to.go <- len - 2*(N.last + 1)
  # The numbers in the next table
  nex <- c(cumsum(last.table)[1:(N.last + 1)], 
                  cumsum(as.numeric(last.table))[(N.last + 2):(N.last + 1 + to.go)] - cumsum(as.numeric(last.table))[1:to.go], 
                  rev(cumsum(as.numeric(last.table))[1:(N.last + 1)]))
  return(nex)
}
