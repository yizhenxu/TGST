#' Cross Validation
#'
#' This function allows you to compute the average of misdiagnoses rate for viral failure and the optimal risk under min-\eqn{\lambda} rules
#' from K-fold cross-validation.
#' @param Obj An object of class TGST. 
#' @param lambda A user-specified weight that reflects relative loss for the two types of misdiagnoses, taking value in \eqn{[0,1]}. \eqn{Loss=\lambda*I(FN)+(1-\lambda)*I(FP)}.
#' @param K Number of folds in cross validation. The default is 10.
#' @return  Cross-validation results.
#' @keywords Cross validation, optimal risk, FNR, FPR.
#' @export
#' @examples
#' d = Simdata
#' Z = d$Z # True Disease Status
#' S = d$S # Risk Score
#' phi = 0.1 #10% of patients taking viral load test
#' Obj = TGST(Z, S, phi, method="nonpar")
#' lambda = 0.8
#' CV.TGST(Obj, lambda, K=10)

CV.TGST <- function(Obj, lambda, K = 10){
  Z <- Obj@Z
  S <- Obj@S
  phi <- Obj@phi
  data <- cbind(Z,S)
  n <- length(Z)
  n.sub <- floor(n/K)
  r <- n-n.sub*K
  if(r==0){
    fld.flg <- sample(rep(1:K,n.sub))
  } else {
    fld.flg <- sample(c(rep(1:K,n.sub),1:r))
  }
  
  fnr.fpr <- matrix(NA,ncol=2,nrow=K)
  risk <- rep(NA,K)
  
  for(i in 1:K){
    
    train.dat <- data[fld.flg!=i,]
    val.dat <- data[fld.flg==i,]
    p <- mean(val.dat[,1])
    
    if(Obj@Nonparametric==TRUE){
      
      opt.rule <- Opt.nonpar.rule(train.dat[,1],train.dat[,2],phi,lambda)
      fnr.fpr[i,] <- nonpar.fnr.fpr(val.dat[,1],val.dat[,2],opt.rule[1],opt.rule[2])
      
    } else {
      
      opt.rule <- Opt.semipar.rule(train.dat[,1],train.dat[,2],phi,lambda)
      fnr.fpr[i,] <- semipar.fnr.fpr(val.dat[,1],val.dat[,2],opt.rule[1],opt.rule[2])
      
    } 
    
    risk[i] <- fnr.fpr[i,1]*p*lambda + fnr.fpr[i,2]*(1-p)*(1-lambda)
    
  }
  
  result <- apply(cbind(fnr.fpr,risk),2,mean)
  names(result) <- c("avr.FNR","avr.FPR","avr.risk")
  
  return(result)
  
}