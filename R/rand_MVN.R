#' Simulate N realizations of X ~ MVN(mu,Sigma)
#'
#' @param mu A n by 1 column vector
#' @param Sigma A n by n matrix
#' @param N An integer
#'
#' @export
rand_MVN <- function(mu, Sigma, N) {
  sol <- c()
  L <- chol(Sigma)
  for(i in 1:N) {
    z <- rnorm(length(mu), mean = 0, sd = 1)
    z <- mu + (L%*%z)
    sol <- cbind(sol,z)
  }
  return(sol)
}
