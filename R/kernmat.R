
kernmat <- function(kf, mat1, mat2)
{
  # Alternative function to \code{\link{kernlab::kernelMatrix}} for two matrices \code{mat1,mat2} with observations in rows
  # @note This function uses double for loops in R and is less cumbersome when "\code{mat2}" is not equal to "\code{mat1}"
  
  if (requireNamespace("kernlab", quietly = TRUE)) {
    res1 <- kernlab::kernelMatrix(kf, mat1, mat2)
  } else {
    stopifnot(is.matrix(mat1) && is.matrix(mat2))
    res1 <- matrix(0, nrow = nrow(mat1), ncol = nrow(mat2))
    dimnames(res1) <- list(rownames(mat1), rownames(mat2))
    for (i in seq(nrow(mat1))) {
      for(j in seq(nrow(mat2))) {
        res1[i,j] <- kf(mat1[i, ], mat2[j, ])
      }
    }
  }
  return(res1)
}
