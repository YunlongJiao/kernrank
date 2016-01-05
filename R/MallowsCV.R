#' Fits a Mallows mixture model to ranking data.
#' 
#' Fits the Mallows mixture model to full ranking data, using
#' Kendall distance and an EM algorithm, for clustering permutations.
#' 
#' 
#' @param datas Matrix of fully-ranked data or permutations.
#' @param G Number of modes, 2 or greater.
#' @param weights Vector of frequencies of each permutation observed,
#' which cannot contain 0 and should be of equal length with nrow of datas since cv split is done by partitioning weights!
#' @param ... Arguments passed to main function Mallows.
#' @param nfolds nfold-fold cv created each time.
#' @param nrepeats cv repeated nrepeats times.
#' @param ntry Number of random initializations to restart each run and best fit with max likelihood chosen.
#' @param logsumexp.trick Logical. Use logsumexp trick to compute log-likelihood.
#' @return See output of FormatOut.
#' @author Yunlong Jiao
#' @importFrom caret createMultiFolds
#' @export
#' @keywords cluster Mallows mixture
#' 

MallowsCV <- function(datas, G, weights = NULL, ..., seed=26921332, nfolds=5, nrepeats=10, ntry = 3, logsumexp.trick = TRUE)
{
  # @param weights 
  # @param nfolds,nrepeats create nfold-fold cv repeated nrepeats times
  # @param ntry algorithm running ntry times for training set trying to avoid local maxima
  
  if (!is.null(seed)) set.seed(seed)
  if (is.null(weights)) weights <- rep(1, nrow(datas))
  
  abils <- ncol(datas)
  perm <- do.call("rbind", permn(abils))
  perm.info <- KendallInfo(perm)
  
  nsample <- sum(weights)
  foldIndices <- caret::createMultiFolds(1:nsample, k=nfolds, times=nrepeats)
  
  wcumsum <- c(0, cumsum(weights))
  weightIndices <- lapply(foldIndices, function(i){
    tab <- table(cut(i, wcumsum, right = TRUE, include.lowest = FALSE))
    as.vector(unname(tab))
  })
  
  foldres <- lapply(weightIndices, function(wtr){
    res <- lapply(seq(ntry), function(i){
      Mallows(datas = datas, G = G, weights = wtr, ..., logsumexp.trick = logsumexp.trick)
    })
    
    # pick up best fit
    likes <- sapply(res, function(model) model[["min.like"]][max(which(model[["min.like"]] != 0))])
    res <- res[[which.max(likes)]]
    
    # calculate cross-validate loglik
    cent <- do.call("rbind", res$R)
    if (!grepl("kernel", res$key)) {
      cent <- KendallInfo(cent)
    }
    dists.to.Rg <- distL2(r = perm.info, centers = cent)
    C.lam <- C_lam(res$lambda, dists.to.Rg = dists.to.Rg, return.logC = logsumexp.trick)
    zcol <- grep("pvals", colnames(res$datas))
    acol <- grep("dists", colnames(res$datas))
    res$cv.loglik <- Likelihood(z = res$datas[ ,zcol], p = res$p, C.lam = C.lam, lambda = res$lambda, 
                                all.dists.data = res$datas[ ,acol], weights = weights-wtr, use.logC = logsumexp.trick)
    res
  })
  
  foldres
}
