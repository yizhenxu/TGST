#' Nonparametric ROC Analysis
#'
#' This function performs ROC analysis for tripartite rules by nonparametric approach. If \eqn{plot=TRUE}, the ROC curve is returned.
#' @param Z True disease status (No disease / treatment success coded as Z=0, diseased / treatment failure coded as Z=1). 
#' @param S Risk score. 
#' @param phi Percentage of patients taking viral load test. 
#' @return 
#' AUC The area under the ROC curve.
#' Misdiagnoses rate for viral failure (i.e., false negative rate, FNR) and otherwise (i.e., false positive rate, FPR). 
#' @keywords Nonparametric, ROC, AUC, FNR, FPR.
#' @export
#' @examples
#' data = Simdata
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
  p <- mean(Z,na.rm=TRUE)
  rules <- Rules.set(Z,S,phi)
  n = length(Z)
  ## AUC
  #Write the kth rule in Rule.set as (i_k,j_k), let j_0=0
  #Hphi(u)=argmin_w {G(u)-G(w)<=phi}
  #For Sj in (j_{k-1},j_k], Hphi(Sj)=i_k
  Hphi <- function(Sj,bounds=rules[,1:2]){
    diff <- Sj-bounds[,2]
    diff1 <- diff[diff<=0][1] #Sj-j_k, where Sj in (j_{k-1},j_k]
    indx <- which(diff==diff1)
    return(bounds[indx,1])
  }
  #calculate AUC from eqn (10) pg1177
  auc <- 0
  for(j in 1:n){
    auc <- auc+sum(Z*(1-Z[j])*((S>Hphi(S[j]))+(S==Hphi(S[j]))/2))
  }
  auc <- auc/(n^2*p*(1-p))
  if(plot==TRUE){  
  #ROC curve
  plot(rules[,3],1-rules[,4],type="l",xlab="FNR",ylab="TPR",main="ROC Curve")
  legend('bottomright',paste("AUC=",round(auc,3),sep=" "),bty ="n",cex=0.8)
  }
  outpt = list(AUC=auc, FNR=rules[,3], FPR=rules[,4])
  invisible(outpt)
}

