#' Solve OLS problem using SVD decomposition
#'
#' @param X A n by d matrix
#' @param y A n by 1 column vector
#'
#' @export
gd_logistic_regression <- function(X,y,alpha,tol) {
  llhood <- c()
  d <- dim(X)
  M <- diag(d[2])
  beta <- 0*diag(M)
  err <- 1
  count <- 0
  while(err > tol && count < 1000) {
    llhood <- c(llhood,t(y) %*% X %*% beta -sum(log(1/(1+exp(X%*%beta)))))
    gradllhood <- t(X)%*%(y - 1/(1+exp(-X%*%beta)))
    beta <- beta - alpha*M%*%gradllhood
    err <- norm(- alpha*M%*%gradllhood)
    count <- count + 1
  }
  return(list(beta,llhood))
}
