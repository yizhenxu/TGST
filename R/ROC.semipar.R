#' Semiparametric ROC Analysis
#'
#' This function performs ROC analysis on the rules from nonparametric approach. If \eqn{plot=TRUE}, the ROC curve is returned.
#' @param Z True disease status (No disease / treatment success coded as Z=0, diseased / treatment failure coded as Z=1). 
#' @param S Risk score. 
#' @param phi Percentage of patients taking viral load test. 
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

semipar.fnr.fpr <- function(Z,S,l,u){
  if(length(l)!=length(u))   #l is lower cutoff, u is upper cutoff
    cat("***** Warning: Wrong rules set. \n")
  p <- mean(Z)
  temp <- density(S) #marginal density (S)
  fit <- glm(Z~ S, family=binomial)
  beta0star <- fit$coef[1]-log(p/(1-p))
  t <- exp(beta0star+temp$x*fit$coef[2]) #g1=t*g0 under exp tilt assumption
  g1 <- temp$y/(p+(1-p)/t)
  g0 <- temp$y/(p*t+1-p)
  x <- temp$x
  
  len <- length(x)
  dif <- x[2:len]-x[1:(len-1)]
  
  cal.fnr <- function(dens,a){  
    if( a>max(x) ){
      area <- 1
    } else if( a<min(x) ){
      area <- 0
    } else {
      diff <- a-x
      diff1 <- diff[diff<=0][1] 
      indx <- which(diff==diff1)#return index of nearest right endpoint
      area <- sum(dens[1:(indx-2)]*dif[1:(indx-2)])+dens[indx-1]*(a-x[indx-1])
    }
    return(area)
  } 
  
  fnr.fpr <- NULL
  K <- length(l)
  for( i in 1:K){
    fnr <- cal.fnr(g1,l[i])
    fpr <- 1-cal.fnr(g0,u[i])
    fnr.fpr <- rbind(fnr.fpr,c(fnr,fpr))        
  }
  
  return(fnr.fpr)

}



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
