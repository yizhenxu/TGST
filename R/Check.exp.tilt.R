#' Check exponential tilt model assumption
#'
#' This function provides graphical assessment to the suitability of the exponential tilt model for risk score in finding optimal tripartite rules by semiparametric approach. 
#' \deqn{g1(s)=exp(\tilde{\beta}_{0}+\beta_{1}*s)*g0(s)} 
#' @param Z True disease status (No disease / treatment success coded as Z=0, diseased / treatment failure coded as Z=1). 
#' @param S Risk score. 
#' @return 
#' Plot of empirical density for risk score S, joint empirical density for (S,Z=1) and (S,Z=0), and the density under the exponential tilt model assumption for (S,Z=1) and (S,Z=0).
#' @keywords Semiparametric, exponential tilt model.
#' @import ggplot2
#' @export
#' @examples
#' data = Simdata
#' Z = d$Z # True Disease Status
#' S = d$S # Risk Score
#' Check.exp.tilt( Z, S)

Check.exp.tilt <- function(Z,S){
  
  ### graphical assessment to exponential tilt model ###
  data <- cbind(Z,S)
  Z <- data[complete.cases(data),1]
  S <- data[complete.cases(data),2]
  p <- mean(Z,na.rm=TRUE)
  fit <- glm(Z~ S, family=binomial)
  
  
  temp <- density(S) #marginal density (S)
  vec <- data.frame(x=temp$x,y=temp$y,Lines=rep("S",length(temp$x)))
  
  temp1 <- density(S[Z==1]) #conditional density (S|Z=1)
  temp1$y <- temp1$y*p #joint density (S,Z=1)
  vec1 <- data.frame(x=temp1$x,y=temp1$y,Lines=rep("S,Z=1 Empirical",length(temp1$x)))
  
  temp0 <- density(S[Z==0]) #conditional density (S|Z=0)
  temp0$y <- temp0$y*(1-p) #joint density (S,Z=0)
  vec0 <- data.frame(x=temp0$x,y=temp0$y,Lines=rep("S,Z=0 Empirical",length(temp0$x)))
  
  beta0star <- fit$coef[1]-log(p/(1-p))
  t <- exp(beta0star+temp$x*fit$coef[2]) #g1=t*g0 under exp tilt assumption
  #by g = p*g1+(1-p)*g0 = [p*t+(1-p)]*g0 = [p+(1-p)/t]*g1
  #plot the joint density (S,Z=1), (S,Z=0) under exponential tilt model assumption
  g1 <- temp$y/(p+(1-p)/t)
  g0 <- temp$y/(p*t+1-p)
  
  vec1e <- data.frame(x=temp$x,y=g1*p,Lines=rep("S,Z=1 Exp Tilt",length(temp$x)))
  vec0e <- data.frame(x=temp$x,y=g0*(1-p),Lines=rep("S,Z=0 Exp Tilt",length(temp$x)))
  
  dat <- rbind(vec,vec1,vec1e,vec0,vec0e)
  colors <- c("red","blue","blue","green","green")
  types <- c("solid","solid", "dashed","solid","dashed")
  library(ggplot2)
  ggplot(dat, aes(x=x, y=y, colour = Lines,group=Lines,linetype = Lines)) + ggtitle("Check Exponential Tilt Model Assumption") + 
    geom_line() +
    scale_colour_manual(values =colors ) +
    scale_linetype_manual(values = types)
}
