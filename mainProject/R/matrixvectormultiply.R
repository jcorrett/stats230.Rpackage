#' Compute matrix-vector product ABx
#'
#' @param A A matrix
#' @param B A matrix
#' @param x A column vector
#' @param switch_var A boolean
#'
#' @export
matrixvectormultiply <- function(A, B, x, switch_var) {
  if (switch_var == TRUE) {
    prod <- A%*%B
    prod <- prod%*%x
  } else {
    prod <- B%*%x
    prod <- A%*%prod
  }
  return(prod)
}
