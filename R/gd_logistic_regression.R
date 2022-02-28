#' Solve logistic regression problem via gradient descent/ascent
#'
#' @param X A n by d matrix
#' @param y A n by 1 column vector
#' @param alpha A positive integer learning rate
#' @param tol A non-negative error tolerance
#'
#' @export
gd_logistic_regression <- function(X,y,alpha,tol) {
  llhood <- c()
  d <- dim(X)
  M <- diag(d[2])
  beta <- as.matrix(rnorm(d[2]))
  err <- Inf #To ensure we enter the loop
  while(err > tol) {
    llhood <- c(llhood,t(y) %*% X %*% beta -sum(log(1+exp(X%*%beta))))
    gradllhood <- t(X)%*%(1/(1+exp(-X%*%beta)) - y)
    beta <- beta - alpha*M%*%gradllhood
    err <- norm(- alpha*M%*%gradllhood)
  }
  p <- 1/(1+exp(-X%*%beta))
  W <- diag(as.vector(p*(1-p)))
  CI_bounds <- 1.960*sqrt(diag(inv(t(X)%*%W%*%X)))
  return(list(beta,llhood,CI_bounds))
}
