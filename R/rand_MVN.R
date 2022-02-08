#' Simulate N realizations of X ~ MVN(mu,Sigma)
#'
#' @param mu A n by 1 column vector
#' @param Sigma A n by n matrix
#' @param N An integer
#'
#' @export
rand_MVN <- function(mu, Sigma, N) {
  z <- rnorm(N, mean = 0, sd = 1)
  L <- chol(Sigma)
  sol <- mu + (L%*%z)
  return(sol)
}
