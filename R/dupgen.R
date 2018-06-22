#' @importFrom combinat permn
#' 

dupgen <- function(x)
{
  # Generates compatible rankings that breaks tied observations in a ranking,
  # where ties are supposed few and NAs are retained in the compatible rankings.
  # 
  # @param x Vector. 
  # If \code{x} is numeric, the rank vector converted from \code{x} indicate that larger values mean being preferred.
  # 
  # @return A matrix with compatible rankings to "\code{x}" in rows.
  
  stopifnot(is.vector(x))
  if (is.null(names(x)))
    names(x) <- paste0("v", seq_along(x))
  x.prot <- rank(x, na.last = "keep", ties.method = "min")
  
  # Deal with unique values or NA
  id.uniq <- is.na(x.prot) | (!duplicated(x.prot, fromLast=FALSE) & !duplicated(x.prot, fromLast=TRUE))
  res <- t(as.matrix(x.prot[id.uniq])) # nrow set to 1
  
  # Deal with duplicates
  x.tab <- table(x.prot)
  f = factor(x.prot, levels = as.integer(names(which(x.tab > 1))))
  x.sp <- split(x.prot, f)
  for (dup in x.sp) {
    t <- combinat::permn((dup[1]):(dup[1]+length(dup)-1))
    t <- do.call("rbind", t)
    colnames(t) <- names(dup)
    res <- cbind(res[rep(seq_len(nrow(res)), times = nrow(t)), , drop=F],
                 t[rep(seq_len(nrow(t)), each = nrow(res)), , drop=F])
  }
  
  return(res[ , names(x)])
}
