#' @importFrom combinat permn
#' @importFrom stats runif
#' 

MallowsKernel <- function(datas, G, iter, weights, eqlam, tol = 1e-3, iterin = iter, key = "kernelMallows", 
                          logsumexp.trick = TRUE, exh = FALSE)
{
  # Feature points for data
  data.info <- KendallInfo(datas)
  # Dimension of feature space
  dfs <- ncol(data.info)
  
  # Number of subjects
  N <- nrow(datas)
  # Number of items being ranked
  abils <- ncol(datas)
  
  # Generate all full permutations
  # DO NOT RUN for over-sized permutation
  if (abils > 8)
    warning("Size of permutations is very large, this might take a while...")
  perm <- do.call("rbind", combinat::permn(abils))
  perm.info <- KendallInfo(perm)
  
  # Initialize the p-value of membership in each cluster
  p <- rep(1/G, G)
  # Initialize fuzzy assignment
  z <- matrix(1/G, nrow = N, ncol = G)
  # Initialize the modal sequences
  R <- lapply(1:G, function(i) runif(dfs, min = -1, max = 1))
  dists.to.Rg <- distL2(r = perm.info, centers = do.call("rbind", R))
  all.dists.data <- distL2(r = data.info, centers = do.call("rbind", R))
  
  # Initialize the lambda values
  lambda <- rep(runif(1, min = 0.1, max = 10), G)
  C.lam <- C_lam(lambda, dists.to.Rg = dists.to.Rg, return.logC = logsumexp.trick)
  
  # Initialize loglike
  like <- Likelihood(z, p, C.lam, lambda, all.dists.data, weights, use.logC = logsumexp.trick)
  
  # cat("Solving...\n")
  likelihood <- 0*(1:iter)
  i <- 1
  
  while (i <= iter) {
    # Save quantities from last iteration
    z1 <- z
    p1 <- p
    R1 <- R
    dists.to.Rg1 <- dists.to.Rg
    all.dists.data1 <- all.dists.data
    lambda1 <- lambda
    C.lam1 <- C.lam
    like1 <- like
    
    # E Step
    z <- EStep(R, NULL, p, lambda, G, N, C.lam, all.dists.data, use.logC = logsumexp.trick)
    
    # M Step
    p <- UpdateP(z, weights)
    
    if (exh) {
      # Centers and lambda are being fully optimized with alternate minimization
      iin <- 1
      while (iin <= iterin) {
        likeold <- like
        
        R <- lapply(1:G, function(i) optR(rold = R[[i]], zz = weights * z[ ,i], lam = lambda[i], data.info = data.info, perm.info = perm.info))
        dists.to.Rg <- distL2(r = perm.info, centers = do.call("rbind", R))
        all.dists.data <- distL2(r = data.info, centers = do.call("rbind", R))
        
        lambda <- UpdateLambda(r = NULL, R = R, z = weights * z, G = G, all.dists.data = all.dists.data, 
                               dists.to.Rg = dists.to.Rg, eqlam = eqlam)
        C.lam <- C_lam(lambda, dists.to.Rg = dists.to.Rg, return.logC = logsumexp.trick)
        
        like <- Likelihood(z, p, C.lam, lambda, all.dists.data, weights, use.logC = logsumexp.trick)
        
        if (is.finite(like) && is.finite(likeold) && abs(like - likeold) < tol) {
          # message("Exh search converged!")
          iin <- iterin
        }
        
        iin <- iin + 1
      }
    } else {
      R <- lapply(1:ncol(z), function(j){
        colSums(data.info * ((weights * z[ ,j]) / sum(weights * z[ ,j])))
      })
      dists.to.Rg <- distL2(r = perm.info, centers = do.call("rbind", R))
      all.dists.data <- distL2(r = data.info, centers = do.call("rbind", R))
      
      lambda <- UpdateLambda(r = NULL, R = R, z = weights * z, G = G, all.dists.data = all.dists.data, 
                             dists.to.Rg = dists.to.Rg, eqlam = eqlam)
      C.lam <- C_lam(lambda, dists.to.Rg = dists.to.Rg, return.logC = logsumexp.trick)
      
      like <- Likelihood(z, p, C.lam, lambda, all.dists.data, weights, use.logC = logsumexp.trick)
    }
    
    likelihood[i] <- like
    
    if (i > 2 && is.finite(likelihood[i]) && is.finite(likelihood[i-1])) {
      if (abs(likelihood[i] - likelihood[i-1]) < tol) {
        # message("Algorithm converged!")
        likelihood[i:iter] <- likelihood[i]
        i <- iter
      } else if (likelihood[i-1] - likelihood[i] > tol) {
        # In case of being stuck at a stationary point corresp to local minima for exhaustive search via alternate minimization
        warning("Jumped out of main loop due to alternating between local minima!")
        z <- z1
        p <- p1
        R <- R1
        dists.to.Rg <- dists.to.Rg1
        all.dists.data <- all.dists.data1
        lambda <- lambda1
        C.lam <- C.lam1
        like <- like1
        likelihood[i:iter] <- likelihood[i-1]
        i <- iter
      }
    }
    
    i <- i + 1
  }
  
  # cat("Formatting output")
  out <- FormatOut(R, p, lambda, z, datas, likelihood, all.dists.data, weights, key)
  return(out)
}
