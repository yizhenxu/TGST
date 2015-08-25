#' Optimal Semiparametric Rule
#'
#' This function gives you the optimal semiparametric tripartite rule that minimizes the min-\eqn{\lambda} rules.
#' @param Z True disease status (No disease / treatment success coded as Z=0, diseased / treatment failure coded as Z=1). 
#' @param S Risk score. 
#' @param phi Percentage of patients taking viral load test. 
#' @param lambda A user-specified weight that reflects relative loss for the two types of misdiagnoses, taking value in \eqn{[0,1]}. \eqn{Loss=\lambda*I(FN)+(1-\lambda)*I(FP)}.
#' @return 
#' Optimal semiparametric rule and its associated misclassification rates (FNR, FPR), optimal lambda risk, and total misclassification rate (TMR). 
#' @keywords Semiparametric, optimal tripartite rules, optimal risk.
#' @export
#' @examples
#' d = Simdata
#' Z = d$Z # True Disease Status
#' S = d$S # Risk Score
#' phi = 0.1 #10% of patients taking viral load test
#' lambda = 0.5
#' Opt.semipar.rule( Z, S, phi, lambda)

Opt.semipar.rule <- function(Z,S,phi,lambda){
  rules <- nonpar.rules(Z,S,phi)
  fnr.fpr <- semipar.fnr.fpr(Z,S,rules[,1],rules[,2])
  fnr <- fnr.fpr[,1]
  fpr <- fnr.fpr[,2]
  p <- mean(Z,na.rm=TRUE)
  risk <- fnr*p*lambda + fpr*(1-p)*(1-lambda)
  index <- which.min(risk)
  opt.risk <- risk[index]
  opt.rule <- rules[index,]
  opt.fnr.fpr <- fnr.fpr[index,]
  TMR <- fnr.fpr[index,1]*p + fnr.fpr[index,2]*(1-p)
  z <- c(opt.rule,opt.fnr.fpr,opt.risk,TMR)
  names(z) <- c("lower.cutoff","upper.cutoff","FNR","FPR","opt.risk","TMR")
  return(z)
}
