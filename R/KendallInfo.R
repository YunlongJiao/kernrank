#' Kendall embedding of pairwise relative ordering
#' 
#' Performs between-column comparison on a matrix of sequences denoting 
#' \code{sign([,i] - [,j])} for \code{i < j}.
#' 
#' 
#' @param r A vector or a matrix of dimension \code{N x n} with sequences in rows.
#' @return A matrix of dimension \code{N x choose(n,2)} with entry values -1/1/0 representing pairwise comparisons of
#' vector values for each row. Specifically, a -1 value denotes that there is 
#' an increase between the two columns, 1 a decrease, and 0 indicates that 
#' the column values are identical in the same row.
#' @note A matrix with one row is returned if the input "\code{r}" is a vector.
#' @author Yunlong Jiao
#' @export
#' @importFrom combinat combn
#' @references
#' Yunlong Jiao, Jean-Philippe Vert. "The Kendall and Mallows Kernels for Permutations." IEEE Transactions on Pattern Analysis and Machine Intelligence (TPAMI), vol. 40, no. 7, pp. 1755-1769, 2018. \href{https://doi.org/10.1109/TPAMI.2017.2719680}{DOI:10.1109/TPAMI.2017.2719680}
#' @keywords Kendall Embedding
#' @examples 
#' r <- do.call('rbind', combinat::permn(1:5))
#' KendallInfo(r)
#' 

KendallInfo <- function(r)
{
  if (is.vector(r)) {
    r <- matrix(r, nrow = 1)
  } else if (!is.matrix(r)) {
    r <- as.matrix(r)
  }
  
  inds <- combinat::combn(ncol(r), 2)
  
  infos <- sign(r[, inds[1, ], drop = FALSE] - r[, inds[2, ], drop =  FALSE])
  return(infos)
}
