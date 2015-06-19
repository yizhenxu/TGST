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
#' data = read.csv("data.csv", header = TRUE) 
#' Z = (data$PVL_Count_Now>1000) #True Status
#' S = data$CD4_count_Now #Risk Score
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
