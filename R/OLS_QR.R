#' Solve OLS problem using QR decomposition
#'
#' @param X A n by d matrix
#' @param y A n by 1 column vector
#'
#' @export
OLS_QR <- function(X,y) {
  QR <- qr(X)
  Q <- qr.Q(QR)
  R <- qr.R(QR)
  beta <- backsolve(R,t(Q)%*%y)
  return(beta)
}
