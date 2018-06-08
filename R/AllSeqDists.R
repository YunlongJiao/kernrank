
AllSeqDists <- function(seqs)
{
  # @note Taken directly from R CRAN package \href{https://cran.r-project.org/web/packages/RMallow/index.html}{RMallow}.
  
  modal <- 1:ncol(seqs)
  infos <- KendallInfo(seqs)
  dists <- apply(infos, 1, function(i) length(which(i == 1)))
  return(dists)
}
