#' Optimal Tripartite Rule
#'
#' \code{OptimalRule} is the main function of \code{TGST} and it gives you the optimal tripartite rule that minimizes the min-\eqn{\lambda} risk based on the type of user selected approach. 
#' The function takes the risk score and true disease status from a training data set and returns the optimal tripartite rule under the specified proportion of patients able to take gold standard test.
#' @param Obj An object of class TGST. 
#' @param lambda A user-specified weight that reflects relative loss for the two types of misdiagnoses, taking value in \eqn{[0,1]}. \eqn{Loss=\lambda*I(FN)+(1-\lambda)*I(FP)}.
#' @return 
#' Optimal tripartite rule. 
#' @keywords Optimal tripartite rules, optimal lambda risk.
#' @import methods
#' @export
#' @examples
#' d = Simdata
#' Z = d$Z # True Disease Status
#' S = d$S # Risk Score
#' phi = 0.1 #10% of patients taking viral load test
#' lambda = 0.5
#' Obj = TGST(Z, S, phi, method="nonpar")
#' OptimalRule(Obj, lambda)


OptimalRule <- function(Obj,lambda){

  p <- mean(Obj@Z,na.rm=TRUE)
  fnr <- Obj@FNR.FPR[,1]
  fpr <- Obj@FNR.FPR[,2]
  risk <- fnr*p*lambda + fpr*(1-p)*(1-lambda)
  index <- which.min(risk)
  if(0){
    opt.risk <- risk[index]
    opt.rule <- Obj@Rules[index,]
    opt.fnr.fpr <- Obj@FNR.FPR[index,]
    TMR <- fnr[index]*p + fpr[index]*(1-p)
    z <- c(opt.rule,opt.fnr.fpr,opt.risk,TMR)
    names(z) <- c("lower.cutoff","upper.cutoff","FNR","FPR","opt.risk","TMR")
    return(z)
  }
  opt.rule <- Obj@Rules[index,]
  names(opt.rule) <- c("lower.cutoff","upper.cutoff")
  result=new("Output", phi=Obj@phi, Z=Obj@Z, S=Obj@S, Rules=Obj@Rules, Nonparametric=Obj@Nonparametric, FNR.FPR=Obj@FNR.FPR, OptRule=opt.rule)
  print(opt.rule)
  return(invisible(result))
  
}
