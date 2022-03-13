#' Metropolis-Hastings Algorithm For simulating E(X^2) for X~Beta(a,b)
#'
#' @param N Number of MCMC steps
#' @param a Beta distribution parameter 1
#' @param b Beta distribution parameter 2
#' @param x0 Initial guess
#' @param sigma Tuning parameter
#'
#' @export
betaX2_MH <- function(N,a,b,x0,sigma) {
  x <- c()
  xcur <- x0
  for(n in 1:N) {
    xprop <- rlnorm(1,meanlog = log(xcur), sdlog = sigma)
    aij <- min(1,((xprop^(2*a-1)*(1-xprop)^(2*b-2))/
                  ((xcur^(2*a-1)*(1-xcur)^(2*b-2)))))
    if(runif(1) <= aij){
      x <- c(x,xprop)
      xcur <- xprop
    } else {
      x <- c(x,xcur)
    }
  }
  return(x)
}
