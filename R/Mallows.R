#' Fits a Mallows mixture model to ranking data
#' 
#' Fits the Mallows mixture model to total rankings, using EM algorithm, 
#' for clustering permutations.
#' 
#' 
#' @param datas Matrix of dimension N x n with sequences in rows.
#' @param G Number of modes, 2 or greater.
#' @param weights Numeric vector of length N denoting frequencies of each permutation observed.
#' Each observation is observed once by default.
#' Notably it must not contain 0 and should be of equal length with nrow(datas).
#' @param iter Maximum number of iterations for EM algorithm.
#' @param iterin Maximum number of iterations for alternate optimization between centers and lambda. 
#' Effective only when performing kernel Mallows with exhaustive optimization.
#' @param tol Stopping precision.
#' @param logsumexp.trick Logical. Whether or not to use log-sum-exp trick to compute log-likelihood.
#' @param seed 
#' @param seed Seed index for reproducible results when optimization is performed. Set to NULL to disable the action.
#' @param key A character string defining the type of Mallows mixture model to perform: 
#' \itemize{
#' \item \emph{copelandMallows} denotes original Mallows mixture model with cluster centers found by Copeland's method
#' \item \emph{bruteMallows} denotes original Mallows mixture model with cluster centers found by brute-force search for optimal Kemeny consensus (not applicable for large n)
#' \item \emph{bordaMallows} denotes original Mallows mixture model with cluster centers as Borda count
#' \item \emph{kernelMallows} denotes kernel version of Mallows mixture model with cluster centers as the barycenter in Euclidean space induced by Kendall embedding
#' \item \emph{kernelGaussian} denotes Gaussian mixture model in the Euclidean space induced by Kendall embedding
#' }
#' @param exhkey A character string. If it greps successfully in "key", an alternate optimization between centers and lambda.
#' Effective only when performing \emph{kernelMallows} with exhaustive optimization.
#' @param eqlamkey A character string. If it greps successfully in "key", the dispersion parameters (lambda) are constrained to be equal for all clusters; otherwise no constraints on lambda.
#' @return List. 
#' \item{key}{Character string indicating the type of Mallows mixture model performed}
#' \item{R}{List of length "G" of cluster centers, each entry being a permutation of length n if original Mallows mixture model is performed, or a numeric vector of length choose(n,2) if kernel version is performed}
#' \item{p}{Numeric vector of length "G" representing the proportion probability of each cluster}
#' \item{lambda}{Numeric vector of length "G" representing the dispersion parameters of each cluster}
#' \item{datas}{A copy of "datas" on which the Mallows mixture model is fitted, combined with "weights", fuzzy assignment membership probability "z", distances to centers in "R"}
#' \item{min.like}{Numeric vector of length "iter" representing fitted likelihood values at each iteration}
#' @author Yunlong Jiao
#' @importFrom combinat permn
#' @export
#' @references 
#' Murphy, T. B., & Martin, D. (2003). Mixtures of distance-based models for ranking data. Computational statistics & data analysis, 41(3), 645-655.
#' @references 
#' Jiao, Y., & Vert, J.-P. (2016). The Kendall and Mallows Kernels for Permutations. 2016. \href{https://hal.archives-ouvertes.fr/hal-01279273}{hal-01279273}
#' @keywords Clustering MallowsMixture
#' @examples 
#' datas <- do.call('rbind', combinat::permn(1:5))
#' G <- 3
#' weights <- runif(nrow(datas))
#' 
#' # fit Mallows mixture model
#' model <- Mallows(datas, G, weights, key = 'bordaMallows')
#' str(model)

Mallows <-function(datas, G, weights = NULL, 
                   iter = 100, iterin = iter, tol = 1e-3, 
                   logsumexp.trick = TRUE, seed = 47631439, 
                   key = c("copelandMallows", "bruteMallows", "bordaMallows", 
                           "kernelMallows", "kernelMallows_Exh", "kernelGaussian", 
                           "copelandMallows_Eqlam", "bruteMallows_Eqlam", "bordaMallows_Eqlam", 
                           "kernelMallows_Eqlam", "kernelMallows_Exh_Eqlam", "kernelGaussian_Eqlam"), 
                   exhkey = "_Exh", eqlamkey = "_Eqlam")
{
  if (!is.null(seed)) 
    set.seed(seed)
  
  if (is.null(weights))
    weights <- rep(1, nrow(datas))
  
  key <- match.arg(key)
  exh <- grepl(exhkey, key)
  eqlam <- grepl(eqlamkey, key)
  if (grepl("kernelMallows", key)) {
    return(MallowsKernel(datas = datas, G = G, iter = iter, weights = weights, eqlam = eqlam, tol = tol, iterin = iterin, exh = exh, logsumexp.trick = logsumexp.trick, key = key))
  } else if (grepl("kernelGaussian", key)) {
    return(GaussianKernel(datas = datas, G = G, iter = iter, weights = weights, eqlam = eqlam, tol = tol, cov.constr = TRUE, key = key))
  }
  
  # Number of subjects.
  N <- nrow(datas)
  # Number of items being ranked.
  abils <- ncol(datas)
  
  if (grepl("bruteMallows", key)) {
    perm <- do.call("rbind", combinat::permn(abils))
  } else {
    perm <- NULL
  }
  
  dists.table <- DistanceDistribution(abils)
  # Initialize the p-value of membership in each cluster.
  p <- rep(1/G, G)
  # Initialize the modal sequences.
  R <- lapply(1:G, function(i) sample(abils))
  # Initialize the lambda values.
  lambda <- rep(runif(1, min = 0.1, max = 10), G)
  C.lam <-  unlist(lapply(lambda, 
                          function(i) C_lam(i, dists.table = dists.table, return.logC = logsumexp.trick)))
  
  # cat("Solving...\n")
  likelihood <- 0*(1:iter)
  infos <- KendallInfo(datas)
  all.dists.data <- AllKendall(r = datas, seqs = do.call("rbind", R), data.info = infos)
  i <- 1
  
  while (i <= iter) {
    # E Step
    z <- EStep(R, datas, p, lambda, G, N, C.lam, all.dists.data, use.logC = logsumexp.trick)
    
    # M Step
    R <- RankAggreg(r = datas, z = weights * z, infos = infos, perm = perm, key = key)
    
    p <- UpdateP(z, weights)
    
    all.dists.data <- AllKendall(r = datas, seqs = do.call("rbind", R), data.info = infos)
    lambda <- UpdateLambda(r = datas, R = R, z = weights * z, G = G, all.dists.data = all.dists.data,
                           dists.table = dists.table, eqlam = eqlam)
    C.lam <- unlist(lapply(lambda, 
                           function(i) C_lam(i, dists.table = dists.table, return.logC = logsumexp.trick)))
    
    likelihood[i] <- Likelihood(z, p, C.lam, lambda, all.dists.data, weights, use.logC = logsumexp.trick)
    if (i > 2 && is.finite(likelihood[i]) && is.finite(likelihood[i-1]) && abs(likelihood[i] - likelihood[i-1]) < tol) {
      # message("Algorithm converged")
      likelihood[i:iter] <- likelihood[i]
      i <- iter
    }
    i <- i + 1
  }
  # cat("Formatting output")
  out <- FormatOut(R, p, lambda, z, datas, likelihood, all.dists.data, weights, key)
  return(out)
}
