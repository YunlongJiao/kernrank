#' @useDynLib kernrank
#' @importFrom Rcpp sourceCpp
#' 

kendall_weight_inner <- function(x, y, key, k, u, normalized, 
                                 eps = 1e-6)
{
  # Core function for the wrap-up function kendall_weight
  
  # Step 1 - Calculate y[order(x)]
  z <- rank(y)[order(x)]
  n <- length(z)
  
  # Step 2 - Calculate non-inversion number by quicksort deduction
  noninv <- kendall_weight_quickR(as.integer(1:n), 
                                  as.integer(z), 
                                  as.integer(n), 
                                  as.integer(key), 
                                  as.integer(k), 
                                  as.numeric(u))
  
  # Step 3 - Output non-inv number
  if (normalized) {
    s <- switch(
      key,
      "1" = n * (n - 1) / 2, # "ken"
      "2" = sum( pmin(1:n,z) * ((n-1):0) ) / n, # "aken"
      "3" = sum( (1:n>=k) & (z>=k) ) * (sum( (1:n>=k) & (z>=k) ) - 1) / 2, # "top"
      "4" = sum(u)^2 + (n - 2) * sum(u * u[z]), # "add"
      "5" = ( sum(u * u[z])^2 - sum(u^2 * u[z]^2) ) / 2 # "mult"
    )
    
    if (abs(s) < eps)
      return(0)
    else
      return(2 * noninv / s - 1)
  } else {
    return(noninv)
  }
}
