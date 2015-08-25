#' Check exponential tilt model assumption
#'
#' This function provides graphical assessment to the suitability of the exponential tilt model for risk score in finding optimal tripartite rules by semiparametric approach. 
#' \deqn{g1(s)=exp(\tilde{\beta0}+\beta1*s)*g0(s)} 
#' @param Z True disease status (No disease / treatment success coded as Z=0, diseased / treatment failure coded as Z=1). 
#' @param S Risk score. 
#' @return 
#' Plot of empirical density for risk score S, joint empirical density for (S,Z=1) and (S,Z=0), and the density under the exponential tilt model assumption for (S,Z=1) and (S,Z=0).
#' @keywords Semiparametric, exponential tilt model.
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
  o<-par("xpd", "mar") #save the original par setting
  layout(rbind(1,2), heights=c(7,1))  # put legend on bottom 1/8th of the chart
  #layout.show(2) #show layout
  # setup margins for the main plot
  par(mar=c(2,4,4,2)+0.1)
  plot(temp, col="red",xlab="",main="Check Exp Tilt Model Assumption")#, xlim=c(0, 100), ylim=c(0, 0.0005))
  temp1 <- density(S[Z==1]) #conditional density (S|Z=1)
  temp1$y <- temp1$y*p #joint density (S,Z=1)
  lines(temp1, col="blue") #plot empirical joint density (S,Z=1)
  temp0 <- density(S[Z==0]) #conditional density (S|Z=0)
  temp0$y <- temp0$y*(1-p) #joint density (S,Z=0)
  lines(temp0, col="green") #plot empirical joint density (S,Z=0)
  beta0star <- fit$coef[1]-log(p/(1-p))
  t <- exp(beta0star+temp$x*fit$coef[2]) #g1=t*g0 under exp tilt assumption
  #by g = p*g1+(1-p)*g0 = [p*t+(1-p)]*g0 = [p+(1-p)/t]*g1
  #plot the joint density (S,Z=1), (S,Z=0) under exponential tilt model assumption
  g1 <- temp$y/(p+(1-p)/t)
  g0 <- temp$y/(p*t+1-p)
  lines(temp$x, g1*p, lty=2, col="blue") #(S,Z=1)
  lines(temp$x, g0*(1-p), lty=2, col="green") #(S,Z=0)
  # setup for no margins on the legend
  par(mar=c(0, 0, 0, 0))
  # c(bottom, left, top, right)
  plot.new()
  legend('center','group',c("S,Z=1 Empirical","S,Z=1 Exp Tilt","S,Z=0 Empirical","S,Z=0 Exp Tilt","S"), lty = c(1,2,1,2,1),
         col=c('blue','blue','green','green','red'),ncol=3,bty ="n",cex=0.8)
  par(o) #back to original par setting
}
