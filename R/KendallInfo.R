#' Kendall embedding of pairwise relative ordering
#' 
#' Performs each column-wise comparison on a matrix of sequences.  A -1 value
#' denotes that there is an increase between the two columns, 1 a decrease, and
#' 0 indicates that the column values are identical in the row.
#' 
#' 
#' @param r Matrix of sequences. Instances in rows.
#' @return Matrix of -1s, 1s, and 0s representing pairwise comparisons of
#' vector values.
#' @author Yunlong Jiao
#' @references Jiao, Yunlong and Vert, Jean-Philippe. The Kendall and Mallows Kernels for Permutations.
#' Proceedings of The 32nd International Conference on Machine Learning (ICML), JMLR: W&CP 37, 1935-1944, 2015.
#' @keywords Kendall Distance

KendallInfo <- function(r)
{
  inds <- combn(ncol(r), 2)
  
  if (is.vector(r)) {
    r <- as.matrix(r, nrow = 1)
  } else if (!is.matrix(r)) {
    r <- as.matrix(r)
    attr(r, "dimnames") <- NULL
  }
  
  infos <- sign(r[, inds[1, ], drop = FALSE] - r[, inds[2, ], drop =  FALSE])
  return(infos)
}
