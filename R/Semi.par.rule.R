#' Semiparametric Rule
#'
#' This function gives you the optimal semiparametric tripartite rule that minimizes TMR (total misclassification risk).  
#' @param Z True disease status (No disease / treatment success coded as Z=0, diseased / treatment failure coded as Z=1). 
#' @param S Risk score. 
#' @param phi Percentage of patients taking viral load test. 
#' @return 
#' Semiparametric rule and its associated TMR.
#' @keywords Semiparametric, optimal tripartite rules, optimal risk.
#' @export
#' @examples
#' data = Simdata
#' Z = d$Z # True Disease Status
#' S = d$S # Risk Score
#' phi = 0.1 #10% of patients taking viral load test
#' Semi.par.rule( Z, S, phi)

Semi.par.rule <- function(Z,S,phi){

    fit <- glm(Z~ S, family=binomial)
    ####
    hosmerlem <- function (y, yhat, g = 10) 
      {
        cutyhat <- cut(yhat, breaks = quantile(yhat, probs = seq(0, 1, 1/g)), include.lowest = T)
        obs <- xtabs(cbind(1 - y, y) ~ cutyhat)
        expect <- xtabs(cbind(1 - yhat, yhat) ~ cutyhat)
        chisq <- sum((obs - expect)^2/expect)
        P <- 1 - pchisq(chisq, g - 2)
        c("X^2" = chisq, Df = g - 2, "P(>Chi)" = P)
    }
    #yhat <- fit$fitted
    #y <- Z[is.na(Z)==FALSE]
    #hosmerlem(y,yhat)
    ####
    center=-fit$coef[1]/fit$coef[2]
    dist.abs <- unique(sort(abs(S-center)))#all possible half widths
    dist.n <- length(dist.abs)
    for(i in 1:dist.n){
      dist.prob <- mean((S<=(center+dist.abs[i])) * (S>=(center-dist.abs[i])),na.rm=TRUE)
      if(dist.prob>phi){
        ind=i-1
        break
      }
    }
    if(phi>0)
      half.width <- dist.abs[ind]
#    dist.prob <- NULL
#    for(i in 1:dist.n){
#      dist.prob <- c(dist.prob, mean((S<=(center+dist.abs[i])) * (S>=(center-dist.abs[i]))))
#    }
#    if(phi>0)
#      half.width <- dist.abs[sum(dist.prob<=phi)]
    if(phi==0) 
      half.width <- 0
    cutoff <- c(center-half.width, center+half.width)
    fnr.fpr <- c(mean((S<cutoff[1])*Z,na.rm=TRUE)/mean(Z,na.rm=TRUE),
                                mean((S>cutoff[2])*(1-Z),na.rm=TRUE)/mean(1-Z,na.rm=TRUE))
    p <- mean(Z,na.rm=TRUE)
    TMR <- p*fnr.fpr[1]+(1-p)*fnr.fpr[2]
    outpt <- c(cutoff,TMR)
    names(outpt) <- c("lower.cutoff", "upper.cutoff", "minTMR")
    return(outpt)
}