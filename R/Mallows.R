#' Fits a Mallows mixture model to ranking data
#' 
#' Fits the Mallows mixture model to total rankings, using EM algorithm, 
#' for clustering permutations.
#' 
#' 
#' @param datas Matrix of dimension \code{N x n} with sequences in rows.
#' @param G Number of modes, 2 or greater.
#' @param weights Numeric vector of length \code{N} denoting frequencies of each permutation observed.
#' Each observation is observed once by default.
#' Notably it must not contain 0 and should be of equal length with \code{nrow(datas)}.
#' @param iter Maximum number of iterations for EM algorithm.
#' @param iterin Maximum number of iterations for alternate optimization between centers and lambda. 
#' Effective only when performing kernel Mallows with exhaustive optimization.
#' @param tol Stopping precision.
#' @param logsumexp.trick Logical. Whether or not to use log-sum-exp trick to compute log-likelihood.
#' @param seed Seed index for reproducible results when optimization is performed. Set to NULL to disable the action.
#' @param key A character string defining the type of Mallows mixture model to perform: 
#' \itemize{
#' \item \code{copelandMallows} denotes original Mallows mixture model with cluster centers found by Copeland's method
#' \item \code{bruteMallows} denotes original Mallows mixture model with cluster centers found by brute-force search for optimal Kemeny consensus (not applicable for large n)
#' \item \code{bordaMallows} denotes original Mallows mixture model with cluster centers as Borda count
#' \item \code{kernelMallows} denotes kernel version of Mallows mixture model with cluster centers as the barycenter in Euclidean space induced by Kendall embedding
#' \item \code{kernelGaussian} denotes Gaussian mixture model in the Euclidean space induced by Kendall embedding
#' }
#' @param exhkey DO NOT CHANGE. A character string. If it greps successfully in "\code{key}", an alternate optimization between centers and lambda.
#' Effective only when performing "\code{kernelMallows}" with exhaustive optimization.
#' @param eqlamkey DO NOT CHANGE. A character string. If it greps successfully in "\code{key}", the dispersion parameters (or lambda) are constrained to be equal for all clusters; otherwise no constraints on lambda.
#' @return List. 
#' \item{key}{Character string indicating the type of Mallows mixture model performed}
#' \item{R}{List of length "\code{G}" of cluster centers, each entry being a permutation of length "\code{n}" if original Mallows mixture model is performed, or a numeric vector of length \code{choose(n,2)} if kernel version is performed}
#' \item{p}{Numeric vector of length "\code{G}" representing the proportion probability of each cluster}
#' \item{lambda}{Numeric vector of length "\code{G}" representing the dispersion parameters of each cluster}
#' \item{datas}{A copy of "\code{datas}" on which the Mallows mixture model is fitted, combined with "\code{weights}", fuzzy assignment membership probability "\code{z}", distances to centers in "\code{R}"}
#' \item{min.like}{Numeric vector of length "\code{iter}" representing fitted likelihood values at each iteration}
#' @author Yunlong Jiao
#' @importFrom combinat permn
#' @importFrom stats runif
#' @export
#' @seealso \code{\link{MallowsCV}}
#' @references 
#' Thomas Brendan Murphy, Donal Martin. "Mixtures of distance-based models for ranking data." Computational Statistics & Data Analysis, vol. 41, no. 3, pp. 645-655, 2003. \href{https://doi.org/10.1016/S0167-9473(02)00165-2}{DOI:10.1016/S0167-9473(02)00165-2}
#' @references 
#' Yunlong Jiao, Jean-Philippe Vert. "The Kendall and Mallows Kernels for Permutations." IEEE Transactions on Pattern Analysis and Machine Intelligence (TPAMI), vol. 40, no. 7, pp. 1755-1769, 2018. \href{https://doi.org/10.1109/TPAMI.2017.2719680}{DOI:10.1109/TPAMI.2017.2719680}
#' @keywords Clustering MallowsMixture
#' @examples 
#' datas <- do.call('rbind', combinat::permn(1:5))
#' G <- 3
#' weights <- runif(nrow(datas))
#' 
#' # Fit Mallows mixture model
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
  
  # Number of subjects
  N <- nrow(datas)
  # Number of items being ranked
  abils <- ncol(datas)
  
  if (grepl("bruteMallows", key)) {
    perm <- do.call("rbind", combinat::permn(abils))
  } else {
    perm <- NULL
  }
  
  dists.table <- DistanceDistribution(abils)
  # Initialize the p-value of membership in each cluster
  p <- rep(1/G, G)
  # Initialize the modal sequences
  R <- lapply(1:G, function(i) sample(abils))
  # Initialize the lambda values
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
      # message("Algorithm converged!")
      likelihood[i:iter] <- likelihood[i]
      i <- iter
    }
    i <- i + 1
  }
  # cat("Formatting output")
  out <- FormatOut(R, p, lambda, z, datas, likelihood, all.dists.data, weights, key)
  return(out)
}
