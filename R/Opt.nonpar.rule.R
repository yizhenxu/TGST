#' Optimal Nonparametric Rule
#'
#' This function gives you the optimal nonparametric tripartite rule and the related optimal risk under min-\eqn{\lambda} rules.  
#' @param Z True disease status (No disease / treatment success coded as Z=0, diseased / treatment failure coded as Z=1). 
#' @param S Risk score. 
#' @param phi Percentage of patients taking viral load test. 
#' @param lambda A user-specified weight that reflects relative loss for the two types of misdiagnoses, taking value in \eqn{[0,1]}. \eqn{Loss=\lambda*I(FN)+(1-\lambda)*I(FP)}.
#' @return 
#' Optimal nonparametric rule and its associated risk.
#' @keywords Nonparametric, optimal tripartite rules, optimal risk.
#' @export
#' @examples
#' data = Simdata
#' Z = d$Z # True Disease Status
#' S = d$S # Risk Score
#' phi = 0.1 #10% of patients taking viral load test
#' lambda = 0.5
#' Opt.nonpar.rule( Z, S, phi, lambda)

Opt.nonpar.rule <- function(Z,S,phi,lambda){
  rules <- Rules.set(Z,S,phi)
  p <- mean(Z,na.rm=TRUE)
  risk <- rules[,3]*p*lambda + rules[,4]*(1-p)*(1-lambda)
  index <- which.min(risk)
  opt.risk <- risk[index]
  opt.rule <- rules[index,1:2]
  z <- c(opt.rule,opt.risk)
  names(z)[3] <- "opt.risk"
  return(z)
}
