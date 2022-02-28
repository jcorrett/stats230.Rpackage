#' HMM backward algorithm
#'
#' @param P A transition prob matrix
#' @param E An emission probability matrix
#' @param v An initial state distribution
#' @param y vector of observed states
#' @param s dim of hidden state space
#'
#' @export
backward_HMM <- function(P,E,v,y,s) {
  barray <- matrix(0,nrow = 2,ncol = length(y))
  b <- 1+integer(s)
  for(t in rev(1:(length(y)-1))){
    b_old <- b
    barray[,t+1] <- b_old
    for(i in 1:s){
      b[i] <- 0
      for(j in 1:s){
        b[i] <- b[i]+P[i,j]*E[j,y[t+1]]*b_old[j]
      }
    }
  }
  b_old <- b
  barray[,1] <- b_old
  p <- 0
  for(j in 1:s){
    p <- p+v[j]*E[j,y[1]]*b_old[j]
  }
  return(list(p,barray))
}
