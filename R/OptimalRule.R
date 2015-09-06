#' Optimal Tripartite Rule
#'
#' \code{OptimalRule} is the main function of \code{TVLT} and it gives you the optimal tripartite rule that minimizes the min-\eqn{\lambda} risk based on the type of user selected approach. 
#' The function takes the risk score and true disease status from a training data set and returns the optimal tripartite rule under the specified proportion of patients able to take VL test.
#' @param Obj An object of class TVLT. 
#' @param lambda A user-specified weight that reflects relative loss for the two types of misdiagnoses, taking value in \eqn{[0,1]}. \eqn{Loss=\lambda*I(FN)+(1-\lambda)*I(FP)}.
#' @return 
#' Optimal tripartite rule and its associated misclassification rates (FNR, FPR), optimal lambda risk, and total misclassification rate (TMR). 
#' @keywords Optimal tripartite rules, optimal lambda risk.
#' @export
#' @examples
#' d = Simdata
#' Z = d$Z # True Disease Status
#' S = d$S # Risk Score
#' phi = 0.1 #10% of patients taking viral load test
#' lambda = 0.5
#' Obj = TVLT(Z, S, phi, method="nonpar")
#' OptimalRule(Obj, lambda)


OptimalRule <- function(Obj,lambda){

  p <- mean(Obj@Z,na.rm=TRUE)
  fnr <- Obj@FNR.FPR[,1]
  fpr <- Obj@FNR.FPR[,2]
  risk <- fnr*p*lambda + fpr*(1-p)*(1-lambda)
  index <- which.min(risk)
  opt.risk <- risk[index]
  opt.rule <- Obj@Rules[index,]
  opt.fnr.fpr <- Obj@FNR.FPR[index,]
  TMR <- fnr[index]*p + fpr[index]*(1-p)
  z <- c(opt.rule,opt.fnr.fpr,opt.risk,TMR)
  names(z) <- c("lower.cutoff","upper.cutoff","FNR","FPR","opt.risk","TMR")
  return(z)
  
}
