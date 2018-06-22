
kendall_weight_inner <- function(x, y, method, k, u, normalized)
{
  # Core function for the wrap-up function kendall_weight
  
  # Convert
  
  # Step 1 - Calculate y[order(x)]
  x <- rank(x)
  y <- rank(y)
  z <- as.integer(y[order(x)])
  n <- length(z)
  
  # Step 2 - Calculate non-inversion number by quicksort deduction
  key <- as.integer(switch(
    method,
    "ken" = 1,
    "aken" = 2,
    "top" = 3,
    "add" = 4,
    "mult" = 5
  ))
  noninv <- kendall_weight_quick(1:n, z, n, key, k, u)
  
  # Step 3 - Output non-inv number
  if (normalized) {
    s <- switch(
      method,
      "ken" = n * (n - 1) / 2,
      "aken" = sum( pmin(1:n,z) * ((n-1):0) ) / n,
      "top" = sum( (1:n>=k) & (z>=k) ) * (sum( (1:n>=k) & (z>=k) ) - 1) / 2,
      "add" = sum(u)^2 + (n - 2) * sum(u * u[z]),
      "mult" = ( sum(u * u[z])^2 - sum(u^2 * u[z]^2) ) / 2
    )
    return(2 * noninv / s - 1)
  } else {
    return(noninv)
  }
}
