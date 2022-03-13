#' Metropolis-Hastings Algorithm For simulating E[X^2] for X~Beta(a,b)
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
betaX2_MH <- function(N,a,b,x0,sigma) {
  x <- c()
  xcur <- x0
  for(n in 1:N) {
    xprop <- rlnrom(1,meanlog = xcur, sdlog = sigma)
    a <- min(1,((xprop^(2*a-1)*(1-xprop)^(2*b-2))/
                  ((xcur^(2*a-1)*(1-xcur)^(2*b-2)))))
    if(runif(1) <= a){
      x <- c(x,xprop)
      xcur <- xprop
    } else {
      x <- c(x,xcur)
    }

  }
  return(x)
}
