#' @importFrom combinat permn
#' 

dupgen <- function(x)
{
  # translate tied observations into compatible full rankings
  # ancillary function for duplicating compatible full rankings for partial rankings with tied observations
  # @return compatible rankings in rows of a matrix if ranking contains duplicated ranks
  
  t <- cumsum(c(1, table(x)))
  xupdate <- x * 0
  for (i in seq(length(t))[-1]) {
    id <- which(x == as.numeric(names(t)[i]))
    xupdate[id] <- t[i - 1] # only relative ranks are reserved in xupdate
  }
  
  xfinite <- xupdate[is.finite(xupdate)]
  u <- unique(xfinite[duplicated(xfinite)])
  
  uid <- lapply(u, function(i) which(xupdate == i))
  res <- matrix(xupdate, nrow = 1)
  
  if (length(u) > 0) {
    for (i in seq(length(u))) {
      nl <- length(uid[[i]])
      s <- factorial(nl)
      subres <- matrix(0, nrow = s * nrow(res), ncol = ncol(res))
      subres[ , uid[[i]]] <- (u[i] + do.call("rbind", combinat::permn(nl)) - 1)[rep(1:s, each = nrow(res)), ]
      subres[ , -uid[[i]]] <- (res[ , -uid[[i]], drop = FALSE])[rep(1:nrow(res), times = s), ]
      res <- subres
    }
  }
  res
}
