
UpdateLambda <- function(r, R, z, G, all.dists.data = NULL, dists.to.Rg = NULL, dists.table = NULL, top.bound = 1000, eqlam)
{
  # upper and lower bound for search lambda AND lower/upper bound for expected loglike
  if (!is.null(dists.table)) {
    flower1 <- C_lam(0, dists.table = dists.table, return.logC = FALSE)*sum(dists.table*as.numeric(names(dists.table)))
    flower1 <- rep(flower1, G)
    
    fupper1 <- C_lam(top.bound, 
                     dists.table = dists.table, return.logC = FALSE)*sum(as.numeric(names(dists.table))*dists.table*exp(-top.bound*as.numeric(names(dists.table))))
    fupper1 <- rep(fupper1, G)
    
    likelow <- C_lam(0, dists.table = dists.table, return.logC = TRUE) * colSums(z)
    likeupp <- C_lam(top.bound, 
                     dists.table = dists.table, return.logC = TRUE) * colSums(z) - top.bound * colSums(z * all.dists.data)
  } else if (!is.null(dists.to.Rg)) {
    flower1 <- colSums(dists.to.Rg)/nrow(dists.to.Rg)
    fupper1 <- rowSums(t(dists.to.Rg) * exp(C_lam(rep(top.bound, G), 
                                                  dists.to.Rg = dists.to.Rg, return.logC = TRUE) - top.bound * t(dists.to.Rg)))
    
    likelow <- C_lam(rep(0, G), dists.to.Rg = dists.to.Rg, return.logC = TRUE) * colSums(z)
    likeupp <- C_lam(rep(top.bound, G), dists.to.Rg = dists.to.Rg, return.logC = TRUE) * colSums(z) - top.bound * colSums(z * all.dists.data)
    
    #     while (any(is.na(fupper1))) {
    #       top.bound <- top.bound/2
    #       warning("top.bound reduced to ", top.bound) 
    #       fupper1 <- rowSums(t(dists.to.Rg) * exp(C_lam(rep(top.bound, G), 
    #                                                     dists.to.Rg = dists.to.Rg, return.logC = TRUE) - top.bound * t(dists.to.Rg)))
    #       likeupp <- C_lam(rep(top.bound, G), dists.to.Rg = dists.to.Rg, return.logC = TRUE) * colSums(z) - top.bound * colSums(z * all.dists.data)
    #     }
  } else {
    stop("Unable to calculate Mallows proba distribution normalization constant!")
  }
  
  # Calculate the RHS of the equation needed to determine lambda values
  # AND correct likelow/likeupp for eqlam
  if (eqlam) {
    rhs <- sum(z*all.dists.data)/sum(z) # of length 1!
    likelow <- sum(likelow) # of length 1!
    likeupp <- sum(likeupp) # of length 1!
  } else {
    rhs <- sapply(1:G, function(i) sum(z[, i]*all.dists.data[, i])/sum(z[, i])) # of length G!
  }
  
  if (eqlam && is.null(dists.table) && !is.null(dists.to.Rg)) {
    # START CASE I
    zcol <- colSums(z)
    flower <- sum(zcol * flower1) / sum(zcol) - rhs
    fupper <- sum(zcol * fupper1) / sum(zcol) - rhs
    # Define function for LHS of the equaltion
    LHS2lam <- function(lambda, rhs, dists.to.Rg, zcol){
      ss <- sapply(1:G, function(i) Lambda(lambda = lambda, rhs = 0, dists = dists.to.Rg[ ,i]))
      sum(zcol * ss) / sum(zcol) - rhs
    }
    if (sign(flower) != sign(fupper)) {
      lambda <- uniroot(LHS2lam, interval = c(0, top.bound), 
                        rhs = rhs, dists.to.Rg = dists.to.Rg, zcol = zcol, 
                        f.lower = flower, f.upper = fupper)$root
    } else {
      if (is.finite(likeupp) && likeupp > likelow) {
        warning("Solution Exceeded top bound, defaulting to top.bound = ", top.bound)
        lambda <- top.bound
      } else {
        warning("Solution Exceeded lower bound, defaulting to 0")
        lambda <- 0
      }
    }
    return(rep(lambda, G)) # RETURN CASE I equal lambda and center being any point
  } else {
    # START CASE II/III/IV
    lambda <- 0*(1:G)
    for (i in 1:G) {
      flower  <- flower1[i] - rhs[i]
      fupper <- fupper1[i] - rhs[i]
      # Find the root of the Lambda function, to 
      # determine lambda.
      if (sign(flower) != sign(fupper)) {
        lambda[i] <- uniroot(Lambda, interval = c(0, top.bound), 
                             rhs = rhs[i], dists = dists.to.Rg[ ,i], dists.table = dists.table, 
                             f.lower = flower, f.upper = fupper)$root
      } else {
        if (is.finite(likeupp[i]) && likeupp[i] > likelow[i]) {
          warning("Solution Exceeded top bound, defaulting to top.bound = ", top.bound)
          lambda[i] <- top.bound
        } else {
          warning("Solution Exceeded lower bound, defaulting to 0")
          lambda[i] <- 0
        }
      }
      if (eqlam) {
        return(rep(lambda[i], G)) # RETURN CASE II equal lambda and center being true permutation
      }
    }
    return(lambda) # RETURN CASE III/IV non-equal lambda and center being true permutation or any point
  }
  
  stop("update lambda return error")
}
