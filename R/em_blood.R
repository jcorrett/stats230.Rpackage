#' Solve logistic regression problem via gradient descent/ascent
#'
#' @param n A 4 by 1 vector of blood type counts. n = (nA,nB,nAB,nO)
#' @param p A 3 by 1 vector of allele probabilities. p = (pA,pB,pO)
#' @param tol A non-negative error tolerance
#'
#' @export
em_blood <- function(n,p,tol) {
err <- Inf
  while(err > tol){
    p_old <- p

    mAA <- n[1]*((p[1]^2)/(p[1]^2+2*p[1]*p[3]))
    mAO <- n[1]*((2*p[1]*p[3])/(p[1]^2+2*p[1]*p[3]))
    mBB <- n[2]*((p[2]^2)/(p[2]^2+2*p[2]*p[3]))
    mBO <- n[2]*((2*p[2]*p[3])/(p[2]^2+2*p[2]*p[3]))
    mAB <- n[3]
    mOO <- n[4]

    p <- c((2*mAA+mAO+mAB)/(2*sum(n)), (2*mBB+mBO+mAB)/(2*sum(n)), (2*mOO+mAO+mBO)/(2*sum(n)))
    err <- norm(p-p_old,type="2")
  }
  p_out <- c(p[1]^2,2*p[1]*p[3],p[2]^2,2*p[2]*p[3],2*p[1]*p[2],p[3]^2)

  return(list(p,p_out))
}
