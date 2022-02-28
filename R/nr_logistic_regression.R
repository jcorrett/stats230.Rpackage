#' Solve logistic regression problem via Newton-Raphson algorithm
#'
#' @param X A n by d matrix
#' @param y A n by 1 column vector
#' @param tol A non-negative error tolerance
#'
#' @export
nr_logistic_regression <- function(X,y,tol) {
  llhood <- c()
  d <- dim(X)
  beta <- integer(d[2])
  p <- 1/(1+exp(-X%*%beta))
  err <- Inf #To ensure we enter the loop
  while(err > tol) {
    llhood <- c(llhood,t(y) %*% X %*% beta -sum(log(1+exp(X%*%beta))))
    gradllhood <- t(X)%*%(1/(1+exp(-X%*%beta)) - y)
    M <- t(X)%*%diag(as.vector(p*(1-p)))%*%X

    beta <- beta - solve(M,gradllhood)
    p <- 1/(1+exp(-X%*%beta))
    err <- norm(-solve(M,gradllhood))
  }
  p <- 1/(1+exp(-X%*%beta))
  W <- diag(as.vector(p*(1-p)))
  CI_bounds <- 1.960*sqrt(diag(inv(t(X)%*%W%*%X)))
  return(list(beta,llhood,CI_bounds))
}
