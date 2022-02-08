#' Solve OLS problem using SVD decomposition
#'
#' @param X A n by d matrix
#' @param y A n by 1 column vector
#'
#' @export
OLS_SVD <- function(X,y) {
  SVD <- svd(X)
  beta <- SVD$v%*%diag(1/SVD$d)%*%t(SVD$u)%*%y
  return(beta)
}
