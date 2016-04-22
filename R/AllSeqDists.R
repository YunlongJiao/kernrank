
AllSeqDists <- function(seqs)
{
  # @note Taken directly from \href{https://cran.r-project.org/web/packages/RMallow/index.html}{CRAN package RMallow}.
  
  modal <- 1:ncol(seqs)
  infos <- KendallInfo(seqs)
  dists <- apply(infos, 1, function(i) length(which(i == 1)))
  return(dists)
}
