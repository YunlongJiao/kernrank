#' @importFrom combinat permn
#'

SeqDistribution <- function(N)
{
  # @note Taken directly from R CRAN package \href{https://cran.r-project.org/web/packages/RMallow/index.html}{RMallow}.
  
  seqs <- do.call("rbind", combinat::permn(N))
  all.dists <- AllSeqDists(seqs)
  return(all.dists)
}
