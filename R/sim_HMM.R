#' Simulate HMM
#'
#' @param P A transition prob matrix
#' @param E An emission probability matrix
#' @param v An initial state distribution
#' @param n length of time series (# of steps)
#'
#' @export
sim_HMM <- function(P,E,V,n) {
  x <- c()
  y <- c()

  for(i in 1:n){
    if(i==1){
      x <- c(x,which(rmultinom(1,1,v)==1))
    } else {
      x <- c(x,which(rmultinom(1,1,P[x[i-1],])==1))
    }
    y <- c(y,which(rmultinom(1,1,E[x[i],])==1))

  }
  return(list(x,y))
}
