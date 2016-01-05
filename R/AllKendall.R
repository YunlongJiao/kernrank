#' All Kendall's distances between two sets of rankings.
#' 
#' Calculates all of the Kendall's distances between two different sets of
#' rankings.
#' 

#' 
#' @param r One set of sequences.
#' @param seqs Another set of sequences.
#' @param data.info Optional argument, a -1/1/0 matrix specifying all of the
#' relevant information to calculate Kendall's difference for "r".  Used for
#' efficiency in "Solve".
#' @return Matrix where output[i, j] represents the distance from sequence "i"
#' in "r" to sequence "j" in "seqs".
#' @author Yunlong Jiao
#' @keywords Kendall distance
#' @export
#' @examples
#' ## Not run:
#' 
#' data1 <- do.call("rbind", list(1:5, 5:1, c(3, 2, 1, 4, 5)))
#' data2 <- do.call("rbind", list(1:5, 5:1))
#' # AllKendall(data1, data2)
#' 
#' ## End(Not run)

AllKendall <- function(r, seqs, data.info = NULL)
{
  if (is.null(data.info)) {
    data.info <- KendallInfo(r)
  }
  seqs.info <- KendallInfo(seqs)
  
  dists <- distL2(r = data.info, centers = seqs.info)
  return(dists)
}
