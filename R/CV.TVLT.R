#' Cross Validation
#'
#' This function allows you to compute the average of misdiagnoses rate for viral failure and the optimal risk under min-\eqn{\lambda} rules
#' from K-fold cross-validation.
#' @param Z True disease status (No disease / treatment success coded as Z=0, diseased / treatment failure coded as Z=1). 
#' @param S Risk score. 
#' @param phi Percentage of patients taking viral load test. 
#' @param lambda A user-specified weight that reflects relative loss for the two types of misdiagnoses, taking value in \eqn{[0,1]}. \eqn{Loss=\lambda*I(FN)+(1-\lambda)*I(FP)}.
#' @param K Number of folds in cross validation. The default is 10.
#' @param method The method to be used. The default is "nonpar", which returns result using nonparametric method. Another possible value is "semipar", which returns result estimated by semiparametric method assuming exponential tilt.
#' @return  Cross-validation results.
#' @keywords Cross validation, optimal risk, FNR, FPR.
#' @export
#' @examples
#' data = Simdata
#' Z = d$Z # True Disease Status
#' S = d$S # Risk Score
#' phi = 0.1 #10% of patients taking viral load test
#' lambda =0.5
#' CV.TVLT(Z, S, phi, K = 10, method = "semipar", lambda)

CV.TVLT <- function(Z, S, phi, K = 10, method = "nonpar", lambda){
  or.data <- cbind(Z,S)
  data <- or.data[complete.cases(or.data),]
  Z <- data[,1]
  S <- data[,2]
  n <- length(Z)
  n.sub <- floor(n/K)
  r <- n-n.sub*K
  if(r==0){
    fld.flg <- sample(rep(1:K,n.sub))
  } else {
    fld.flg <- sample(c(rep(1:K,n.sub),1:r))
  }
  
  fpr <- rep(NA,K)
  fnr <- rep(NA,K)
  risk <- rep(NA,K)
  
  for(i in 1:K){
    
    train.dat <- data[fld.flg!=i,]
    val.dat <- data[fld.flg==i,]
    p <- mean(val.dat[,1])
    
    if(method == "nonpar"){
    
      opt.rule <- Opt.nonpar.rule(train.dat[,1],train.dat[,2],phi,lambda)
      fnr[i] <- mean((val.dat[,2]<opt.rule[1])*val.dat[,1])/mean(val.dat[,1])
      fpr[i] <- mean((val.dat[,2]>opt.rule[2])*(1-val.dat[,1]))/mean(1-val.dat[,1])
      risk[i] <- fnr[i]*p*lambda + fpr[i]*(1-p)*(1-lambda)
      
    } else if (method == "semipar"){
      
      opt.rule <- Opt.semipar.rule(train.dat[,1],train.dat[,2],phi,lambda)
  
      temp <- density(val.dat[,2]) #marginal density (S)
      fit <- glm(val.dat[,1]~ val.dat[,2], family=binomial)
      beta0star <- fit$coef[1]-log(p/(1-p))
      t <- exp(beta0star+temp$x*fit$coef[2]) #g1=t*g0 under exp tilt assumption
      g1 <- temp$y/(p+(1-p)/t)
      g0 <- temp$y/(p*t+1-p)
      l <- length(temp$x)
      dif <- temp$x[2:l]-temp$x[1:(l-1)]      
      
      cal.fnr <- function(dens,a){  
        diff <- a-temp$x
        diff1 <- diff[diff<=0][1] 
        indx <- which(diff==diff1)#return index of nearest right endpoint
        area <- sum(dens[1:(indx-2)]*dif[1:(indx-2)])+dens[indx-1]*(a-temp$x[indx-1])
        return(area)
      } 
      
      fnr[i] <- cal.fnr(g1,opt.rule[1])
      fpr[i] <- 1-cal.fnr(g0,opt.rule[2])
      risk[i] <- fnr[i]*p*lambda + fpr[i]*(1-p)*(1-lambda)
      
    } else {
      cat("Wrong method input!\n")
    }
  }
  
  result <- apply(cbind(fnr,fpr,risk),2,mean)
  return(result)
}