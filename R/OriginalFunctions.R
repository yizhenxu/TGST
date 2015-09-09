#' Summary of Disease Status and Risk Score
#'
#' This function allows you to compute the percentage of diseased (equivalent to treatment failure Z=1) and show distribution summary of risk score (S) by the true disease status (Z).
#' @param Z True disease status (No disease / treatment success coded as Z=0, diseased / treatment failure coded as Z=1). 
#' @param S Risk score. 
#' @return 
#' Percentage of treatment failure; 
#' Summary statistics (mean, standard deviation, minimum, median, maximum and IQR) of risk score by true disease status; 
#' Distribution plot.
#' @keywords Prevalence of disease (treatment failure), risk score summary, risk score distribution.
#' @export
#' @examples
#' d = Simdata
#' Z = d$Z # True Disease Status
#' S = d$S # Risk Score
#' TVLT.summ(Z,S)


TVLT.summ = function(Z,S){
  S0 = S[Z==0]
  S1 = S[Z==1]
  #prevalence
  percF = mean(Z,na.rm=TRUE)
  #summary quantile
  summ.S0 = c(mean(S0,na.rm = T),sd(S0,na.rm = T),min(S0,na.rm = T),median(S0,na.rm = T),max(S0,na.rm = T),quantile(S0,probs=c(0.25,0.75),na.rm = T))
  names(summ.S0) = c("mean","sd","min","median","max","IQR.L","IQR.U")
  summ.S1 = c(mean(S1,na.rm = T),sd(S1,na.rm = T),min(S1,na.rm = T),median(S1,na.rm = T),max(S1,na.rm = T),quantile(S1,probs=c(0.25,0.75),na.rm = T))
  names(summ.S1) = c("mean","sd","min","median","max","IQR.L","IQR.U")
  #risk score distribution plot by disease status
  s.min = min(S)
  s.max = max(S)
  dens0 = density(S0, na.rm = T, from = s.min, to = s.max)
  dens1 = density(S1, na.rm = T, from = s.min, to = s.max)
  yrange = range(c(dens0$y,dens1$y))
  plot(dens0,ylim = yrange, col="blue",
       xlab="Risk Score S", ylab="Population Distributions",main="Risk Score Distribution")
  lines(dens1,col="red")
  legend((s.max+s.min)/2,yrange[2], c("Z=0","Z=1"),
         lty=c(1,1),  lwd=c(2.5,2.5),col=c("blue","red")) 
  #output
  z = list(percF=percF,SummaryS0=summ.S0,SummaryS1=summ.S1)
  return(z)
}





#' Optimal Nonparametric Rule
#'
#' This function gives you the optimal nonparametric tripartite rule that minimizes the min-\eqn{\lambda} rules.
#' @param Z True disease status (No disease / treatment success coded as Z=0, diseased / treatment failure coded as Z=1). 
#' @param S Risk score. 
#' @param phi Percentage of patients taking viral load test. 
#' @param lambda A user-specified weight that reflects relative loss for the two types of misdiagnoses, taking value in \eqn{[0,1]}. \eqn{Loss=\lambda*I(FN)+(1-\lambda)*I(FP)}.
#' @return 
#' Optimal nonparametric rule and its associated misclassification rates (FNR, FPR), optimal lambda risk, and total misclassification rate (TMR). 
#' @keywords Nonparametric, optimal tripartite rules, optimal risk.
#' @export
#' @examples
#' d = Simdata
#' Z = d$Z # True Disease Status
#' S = d$S # Risk Score
#' phi = 0.1 #10% of patients taking viral load test
#' lambda = 0.5
#' Opt.nonpar.rule( Z, S, phi, lambda)

Opt.nonpar.rule <- function(Z,S,phi,lambda){
  rules <- nonpar.rules(Z,S,phi)
  fnr.fpr <- nonpar.fnr.fpr(Z,S,rules[,1],rules[,2])
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



#' Nonparametric ROC Analysis
#'
#' This function performs ROC analysis for tripartite rules by nonparametric approach. If \eqn{plot=TRUE}, the ROC curve is returned.
#' @param Z True disease status (No disease / treatment success coded as Z=0, diseased / treatment failure coded as Z=1). 
#' @param S Risk score. 
#' @param phi Percentage of patients taking viral load test. 
#' @param plot Logical parameter indicating if ROC curve should be plotted. Default is \code{plot=TRUE}. If false, then only AUC is calculated.
#' @return 
#' AUC The area under the ROC curve.
#' FNR Misdiagnoses rate for viral failure (false negative rate).
#' FPR Misdiagnoses rate for treatment failure (false positive rate).
#' @keywords Nonparametric, ROC, AUC, FNR, FPR.
#' @export
#' @examples
#' d = Simdata
#' Z = d$Z # True Disease Status
#' S = d$S # Risk Score
#' phi = 0.1 #10% of patients taking viral load test
#' a = ROC.nonpar( Z, S, phi,plot=TRUE)
#' a$AUC
#' a$FNR
#' a$FPR


ROC.nonpar <- function(Z,S,phi,plot=TRUE){
  #consider only complete cases
  data <- cbind(Z,S)
  Z <- data[complete.cases(data),1]
  S <- data[complete.cases(data),2]
  rules <- Rules.set(Z,S,phi)
  
  auc = cal.AUC(Z,S,rules[,1],rules[,2])
  if(plot==TRUE){  
    #ROC curve
    plot(rules[,4],1-rules[,3],type="l",xlab="FPR",ylab="TPR",main="ROC Curve")
    legend('bottomright',paste("AUC=",round(auc,3),sep=" "),bty ="n",cex=0.8)
  }
  outpt = list(AUC=auc, FNR=rules[,3], FPR=rules[,4])
  invisible(outpt)
}





#' Semiparametric ROC Analysis
#'
#' This function performs ROC analysis on the rules from nonparametric approach. If \eqn{plot=TRUE}, the ROC curve is returned.
#' @param Z True disease status (No disease / treatment success coded as Z=0, diseased / treatment failure coded as Z=1). 
#' @param S Risk score. 
#' @param phi Percentage of patients taking viral load test. 
#' @param plot Logical parameter indicating if ROC curve should be plotted. Default is \code{plot=TRUE}. If false, then only AUC is calculated.
#' @return 
#' AUC The area under the ROC curve.
#' FNR Misdiagnoses rate for viral failure (false negative rate).
#' FPR Misdiagnoses rate for treatment failure (false positive rate).
#' @keywords Nonparametric, ROC, AUC, FNR, FPR.
#' @export
#' @examples
#' d = Simdata
#' Z = d$Z # True Disease Status
#' S = d$S # Risk Score
#' phi = 0.1 #10% of patients taking viral load test
#' a = ROC.semipar( Z, S, phi,plot=TRUE)
#' a$AUC
#' a$FNR
#' a$FPR

ROC.semipar <- function(Z,S,phi,plot=TRUE){
  #consider only complete cases
  data <- cbind(Z,S)
  Z <- data[complete.cases(data),1]
  S <- data[complete.cases(data),2]
  
  rules <- nonpar.rules(Z,S,phi)
  
  fnr.fpr <- semipar.fnr.fpr(Z,S,rules[,1],rules[,2])
  
  auc <- cal.AUC(Z,S,rules[,1],rules[,2])
  
  if(plot==TRUE){  
    #ROC curve
    plot(fnr.fpr[,2],1-fnr.fpr[,1],type="l",xlab="FPR",ylab="TPR",main="ROC Curve")
    legend('bottomright',paste("AUC=",round(auc,3),sep=" "),bty ="n",cex=0.8)
  }
  outpt = list(AUC=auc, FNR=fnr.fpr[,1], FPR=fnr.fpr[,2])
  invisible(outpt)  
}
