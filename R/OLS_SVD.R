#' Solve OLS problem using SVD decomposition
#'
#' @param X A n by d matrix
#' @param y A n by 1 column vector
#'
#' @export
OLS_SVD <- function(X,y) {
  SVD <- svd(X)
  U <- SVD$u
  S <- diag(1/SVD$d)
  V <- SVD$v
  beta <- V%*%S%*%t(U)%*%y
  return(beta)
}
