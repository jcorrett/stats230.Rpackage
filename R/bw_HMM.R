#' Baum-Welch Algorithm for HMM
#'
#' @param P A transition prob matrix initial guess
#' @param E An emission probability matrix initial guess
#' @param v An initial state distribution initial guess
#' @param y vector of observed states
#' @param s dim of hidden state space
#' @param ymax dim of observed state space
#' @param tol non-negative convergence criteria
#'
#' @export
bw_HMM <- function(P,E,v,y,s,ymax,tol) {
  a <- forward_HMM(P,E,v,y,s)
  py <- a[[1]]
  a <- a[[2]]
  b <- backward_HMM(P,E,v,y,s)
  b <- b[[2]]
  gamma <- a*b/colSums(a*b,2)
  err <- Inf
  while(err > tol){
    py_old <- py
    P_old <- P
    E_old <- E
    v_old <- v

    v <- gamma[,1]

    P <- 0*P_old
    for(t in 2:length(y)){
      g <- matrix(0,nrow = s,ncol = s)
      for(i in 1:s){
        for(j in 1:s){
          g[i,j] <- b[j,t]*E_old[j,y[t]]*P_old[i,j]*a[i,t-1]
        }
      }
      P <- P + g/sum(g)
    }
    P <- P/rowSums(P)

    for(l in 1:ymax){
      for(i in 1:s){
        Esum <- 0
        for(t in 1:length(y)){
          if(y[t]==l){
            Esum <- Esum + gamma[i,t]
          }
        }
        E[i,l] <- Esum/sum(gamma[i,])
      }
    }

    a <- forward_HMM(P,E,v,y,s)
    py <- a[[1]]
    a <- a[[2]]
    b <- backward_HMM(P,E,v,y,s)
    b <- b[[2]]
    gamma <- a*b/colSums(a*b,2)

    err <- abs(log(py)-log(py_old))
  }


  return(P,E,v)
}
