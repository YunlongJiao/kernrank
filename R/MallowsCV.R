#' Compute cross-validated likelihood for Mallows mixture models
#' 
#' Assess model performance by cross-validated (CV) Mallows likelihood.
#' Do NOT run for large number of ranked alternatives "\code{n}".
#' 
#' 
#' @param datas Matrix of dimension \code{N x n} with sequences in rows.
#' @param G Number of modes, 2 or greater.
#' @param weights Integer vector of length \code{N} denoting frequencies of each permutation observed.
#' Each observation is observed once by default.
#' Notably it must not contain 0 and should be of equal length with \code{nrow(datas)}.
#' @param ... Arguments passed to \code{\link{Mallows}}.
#' @param seed Seed index for reproducible results when creating splits of data for CV.
#' Set to NULL to disable the action.
#' @param nfolds \code{nfold}-fold CV created each time.
#' @param nrepeats CV repeated \code{nrepeats} times.
#' @param ntry Number of random initializations to restart for each CV run. 
#' The best fit returning max likelihood is reported.
#' @param logsumexp.trick Logical. Whether or not to use log-sum-exp trick to compute log-likelihood.
#' @return List of length \code{nfolds x nrepeats}, each entry being the result on each fold containing:
#' \item{...}{See output of \code{\link{Mallows}}}
#' \item{cv.loglik}{Likelihood value assessed against test fold while the mixture model is trained on the training fold}
#' @author Yunlong Jiao
#' @note CV split is done by partitioning "\code{weights}" so that "\code{weights}" must be integers.
#' @importFrom combinat permn
#' @export
#' @seealso \code{\link{Mallows}}
#' @references 
#' Thomas Brendan Murphy, Donal Martin. "Mixtures of distance-based models for ranking data." Computational Statistics & Data Analysis, vol. 41, no. 3, pp. 645-655, 2003. \href{https://doi.org/10.1016/S0167-9473(02)00165-2}{DOI:10.1016/S0167-9473(02)00165-2}
#' @references 
#' Yunlong Jiao, Jean-Philippe Vert. "The Kendall and Mallows Kernels for Permutations." IEEE Transactions on Pattern Analysis and Machine Intelligence (TPAMI), vol. 40, no. 7, pp. 1755-1769, 2018. \href{https://doi.org/10.1109/TPAMI.2017.2719680}{DOI:10.1109/TPAMI.2017.2719680}
#' @keywords Clustering MallowsMixture
#' @examples 
#' datas <- do.call('rbind', combinat::permn(1:5))
#' G <- 3
#' weights <- rbinom(nrow(datas), 100, 0.5) # positive integers
#' 
#' # Cross validate Mallows mixture model
#' cv.model <- MallowsCV(datas, G, weights, key = 'bordaMallows', nfolds = 3, nrepeats = 1)
#' # Averaged cv.loglik over all CV folds
#' mean(sapply(cv.model, function(model) model$cv.loglik))
#' 

MallowsCV <- function(datas, G, weights = NULL, ..., seed = 26921332, nfolds = 5, nrepeats = 10, ntry = 3, logsumexp.trick = TRUE)
{
  # @param weights 
  # @param nfolds,nrepeats create 
  # @param ntry Algorithm running \code{ntry} times for training set trying to avoid local maxima
  
  if (!is.null(seed)) 
    set.seed(seed)
  
  if (is.null(weights)) 
    weights <- rep(1, nrow(datas))
  nsample <- sum(weights)
  
  if (!is.integer(weights) || any(weights <= 0)) 
    stop("Weights must take integers for CV runs!")
  
  abils <- ncol(datas)
  perm <- do.call("rbind", combinat::permn(abils))
  perm.info <- KendallInfo(perm)
  
  if (requireNamespace("caret", quietly = TRUE)) {
    foldIndices <- caret::createMultiFolds(1:nsample, k=nfolds, times=nrepeats)
  } else {
    # This serves as a rather simple alternative to \code{\link{caret::createMultiFolds}}
    foldIndices <- lapply(1:nrepeats, function(i){
      test.fold <- split(sample(nsample, replace = FALSE), rep(1:nfolds, length = nsample))
      train.fold <- lapply(test.fold, function(testid){
        setdiff(1:nsample, testid)
      })
      return(train.fold)
    })
    foldIndices <- unlist(foldIndices, recursive = FALSE, use.names = FALSE)
    names(foldIndices) <- as.vector(outer(
      paste0("Fold", formatC(1:nfolds, width = nchar(nfolds), format = "d", flag = "0")), 
      paste0("Rep", formatC(1:nrepeats, width = nchar(nrepeats), format = "d", flag = "0")), 
      FUN = paste, sep = "."))
  }
  
  wcumsum <- c(0, cumsum(weights))
  weightIndices <- lapply(foldIndices, function(i){
    tab <- table(cut(i, wcumsum, right = TRUE, include.lowest = FALSE))
    as.vector(unname(tab))
  })
  
  foldres <- lapply(weightIndices, function(wtr){
    res <- lapply(seq(ntry), function(i){
      Mallows(datas = datas, G = G, weights = wtr, ..., logsumexp.trick = logsumexp.trick)
    })
    
    # Pick up best fit
    likes <- sapply(res, function(model) model[["min.like"]][max(which(model[["min.like"]] != 0))])
    res <- res[[which.max(likes)]]
    
    # Calculate cross-validate loglik
    cent <- do.call("rbind", res$R)
    if (!grepl("kernel", res$key)) {
      cent <- KendallInfo(cent)
    }
    dists.to.Rg <- distL2(r = perm.info, centers = cent)
    C.lam <- C_lam(res$lambda, dists.to.Rg = dists.to.Rg, return.logC = logsumexp.trick)
    zcol <- grep("pvals", colnames(res$datas))
    acol <- grep("dists", colnames(res$datas))
    res$cv.loglik <- Likelihood(z = res$datas[ ,zcol], p = res$p, C.lam = C.lam, 
                                lambda = res$lambda, all.dists.data = res$datas[ ,acol], 
                                weights = weights - wtr, use.logC = logsumexp.trick)
    res
  })
  
  foldres
}
