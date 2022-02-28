#' Solve logistic regression problem via gradient descent/ascent
#'
#' @param X A n by d matrix
#' @param y A n by 1 column vector
#' @param alpha A positive integer learning rate
#' @param alpha A non-negative error tolerance
#'
#' @export
gd_logistic_regression <- function(X,y,alpha,tol) {
  llhood <- c()
  errors <- c()
  d <- dim(X)
  M <- diag(d[2])
  beta <- as.matrix(rnorm(10))
  err <- 1
  count <- 0
  while(err > tol && count < 100000) {
    llhood <- c(llhood,t(y) %*% X %*% beta -sum(log(1/(1+exp(X%*%beta)))))
    gradllhood <- t(X)%*%(1/(1+exp(-X%*%beta)) - y)
    beta <- beta - alpha*M%*%gradllhood
    err <- norm(- alpha*M%*%gradllhood)
    errors <- c(errors,err)
    count <- count + 1
  }
  return(list(beta,llhood,errors))
}
