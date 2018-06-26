context("Weighted Kendall kernel")
library(kernrank)

n.rep <- 20
n.min <- 5
n.max <- 10

# Errors ------------------------------------------------------------------

n <- sample(n.min:n.max, 1)
x <- rnorm(n); x[1] <- NA
expect_error(kendall_weight(x, x))

n <- sample(n.min:n.max, 1)
x <- rnorm(n); x <- as.matrix(x)
expect_error(kendall_weight(x, x))

# Quick vs naive ----------------------------------------------------------

kendall.naive <- function(x, y){
  stopifnot(!any(is.na(x)) && !any(is.na(y)))
  stopifnot(anyDuplicated(x) == 0 && anyDuplicated(y) == 0)
  x <- rank(x)
  y <- rank(y)

  n <- length(x)
  s <- 0
  for (i in 1:n) {
    for (j in 1:n) {
      if (x[i] < x[j] && y[i] < y[j])
        s <- s + switch(key,
                        "ken" = 1,
                        "aken" = min(x[i],y[i])/n,
                        "top" = ifelse(x[i] >= k && y[i] >= k, 1, 0),
                        "add" = (u[x[i]] + u[x[j]]) * (u[y[i]] + u[y[j]]),
                        "mult" = (u[x[i]] * u[x[j]]) * (u[y[i]] * u[y[j]])
        )
    }
  }
  s
}

kendall.quick <- function(x, y){
  # @note Same pseudo code as kendall_weight but implemented in R instead of C++

  stopifnot(!any(is.na(x)) && !any(is.na(y)))
  stopifnot(anyDuplicated(x) == 0 && anyDuplicated(y) == 0)
  x <- rank(x)
  y <- rank(y)

  # step 1 - y[order(x)]
  z <- y[order(x)]
  n <- length(z)

  # step 2 - directly calculate non-inversion number by quick sort
  noninv <- 0
  func <- function(idx){
    # take also global variables - static z, static u, nonstatic noninv
    # s will be updated after deduction run
    idxhigh <- c() # index associated to values greater than pivot value
    idxlow <- c() # index associated to values less than or equal to pivot value
    # suppose array z is a permutation which has no equal values
    if (length(idx) > 1){
      # randomly pick an index as pivot
      pivot <- sample(idx, 1)

      # update noninv
      cnum <- 0
      cmin <- 0
      ctop <- 0
      cweight <- 0
      cconvweight <- 0
      cmultweight <- 0
      for (i in idx) {
        if (z[i] < z[pivot]) {
          idxlow <- c(idxlow, i)
          cnum <- cnum + 1
          cmin <- cmin + min(i, z[i])/n
          ctop <- ctop + ifelse(i >= k && z[i] >= k, 1, 0)
          cweight <- cweight + u[i]
          cconvweight <- cconvweight + u[z[i]]
          cmultweight <- cmultweight + (u[i]*u[z[i]])
        } else {
          idxhigh <- c(idxhigh, i)
          noninv <<- noninv +
            switch(key,
                   "ken" = cnum,
                   "aken" = cmin,
                   "top" = ifelse(i >= k && z[i] >= k, ctop, 0),
                   "add" = cmultweight + cweight * u[z[i]] + cconvweight * u[i] + cnum * (u[i]*u[z[i]]),
                   "mult" = cmultweight * (u[i]*u[z[i]])
            )
        }
      }

      # deduction
      func(idxhigh)
      func(idxlow)
    }
  }
  # run deduction to get noninv number
  func(1:n)

  # step 3 - output noninv directly
  noninv
}

for (i in seq(n.rep)) {
  n <- sample(n.min:n.max, 1)
  x <- rnorm(n)
  y <- rnorm(n)
  u <- abs(rnorm(n))
  k <- sample(1:n, 1)
  expect_equal(kendall_weight(x, y, method = "top", k = 1, normalized = TRUE),
               kendall_weight(x, y, method = "ken", normalized = TRUE))
  for (key in c("ken", "aken", "top", "add", "mult")) {
    expect_equal(kendall.naive(x, y),
                 kendall.quick(x, y))
    expect_equal(kendall.naive(x, y),
                 kendall_weight(x, y, method = key, k = k, u = u, normalized = FALSE))
    
    ss <- kendall_weight(x, x, method = key, k = k, u = u, normalized = TRUE)
    if (abs(ss) >= 1e-6) expect_equal(ss, 1)
    
    ss <- kendall_weight(x, order(rev(order(x))), method = key, k = k, u = u, normalized = TRUE)
    if (abs(ss) >= 1e-6) expect_equal(ss, -1)
  }
}
