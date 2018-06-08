
kendall_partial_inner <- function(x,y)
{
  # core function for the wrap-up function kendall_partial
  
  if(length(x) != length(y)){
    stop('Vectors have different lengths!')
  } else {
    n <- length(x)
  }
  
  i <- which(!is.na(x)); k <- length(i)
  j <- which(!is.na(y)); m <- length(j)
  if(k==1 || m==1)
    stop('Length of observed entries should be at least 2!')
  
  # find 4 cases
  idxcommon <- intersect(i,j); l0 <- length(idxcommon)
  idxonlyx <- setdiff(i,j); l1 <- length(idxonlyx)
  idxonlyy <- setdiff(j,i); l2 <- length(idxonlyy)
  idxcomp <- setdiff(1:n,union(i,j)); l3 <- length(idxcomp)
  
  icommon <- which(i %in% idxcommon, arr.ind = TRUE); ionlyx <- setdiff(1:k,icommon)
  jcommon <- which(j %in% idxcommon, arr.ind = TRUE); jonlyy <- setdiff(1:m,jcommon)
  
  if(l0 > 0 && (l1 > 0 || l2 > 0)){
    rxcommon <- rank(x[idxcommon])
    rycommon <- rank(y[idxcommon])
  }
  
  if((l1 > 0 && l2 > 0) || (l0 > 0 && (l1 > 0 || l2 > 0 || l3 > 0))){
    rx <- rank(x[i])
    ry <- rank(y[j])
  }
  
  # 5 cases
  if(l0 < 2){
    ka <- 0
  } else {
    # cat('case a holds!\n')
    ka <- kendall_corr(x[idxcommon],y[idxcommon])*choose(l0,2)/choose(n,2)
  }
  
  if(l0 == 0 || l1 == 0){
    kc <- 0
  } else {
    # cat('case c holds!\n')
    kc <- sum((2*(rx[icommon]-rxcommon)-k+l0)*(2*ry[jcommon]-m-1))/choose(n,2)/(m+1)
  }
  
  if(l0 == 0 || l2 == 0){
    kd <- 0
  } else {
    # cat('case d holds!\n')
    kd <- sum((2*(ry[jcommon]-rycommon)-m+l0)*(2*rx[icommon]-k-1))/choose(n,2)/(k+1)
  }
  
  if(l0 == 0 || l3 == 0){
    ke <- 0
  } else {
    # cat('case e holds!\n')
    ke <- l3*sum((2*rx[icommon]-k-1)*(2*ry[jcommon]-m-1))/(k+1)/(m+1)/choose(n,2)
  }
  
  if(l1 == 0 || l2 == 0){
    kg <- 0
  } else {
    # cat('case g holds!\n')
    kg <- (-1)*sum(2*rx[ionlyx]-k-1)*sum(2*ry[jonlyy]-m-1)/(k+1)/(m+1)/choose(n,2)
  }
  
  return(ka+kc+kd+ke+kg)
}
