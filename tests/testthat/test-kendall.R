context("Kendall kernel")
library(kernrank)

n.rep <- 100
n.min <- 5
n.max <- 20

# Errors ------------------------------------------------------------------

n <- sample(n.min:n.max, 1)
x <- rnorm(n); x[1] <- NA
expect_error(kendall_total(x, x))

n <- sample(n.min:n.max, 1)
x <- rnorm(n); x <- as.matrix(x)
expect_error(kendall_total(x, x))

# cor vs cor.fk -----------------------------------------------------------

for (i in seq(n.rep)) {
  n <- sample(n.min:n.max, 1)
  x <- rnorm(n)
  y <- rnorm(n)
  expect_equal(cor(x, y, method = "kendall"), kendall_total(x, y))
  expect_equal(cor(x, y, method = "kendall"), pcaPP::cor.fk(x, y))
}

for (i in seq(n.rep)) {
  n <- sample(n.min:n.max, 1)
  x <- rep(1:5, length.out=n)[sample(n)]
  y <- rep(1:5, length.out=n)[sample(n)]
  expect_equal(cor(x, y, method = "kendall"), kendall_total(x, y))
  expect_equal(cor(x, y, method = "kendall"), pcaPP::cor.fk(x, y))
}
