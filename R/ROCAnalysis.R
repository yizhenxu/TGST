#' ROC Analysis
#'
#' This function performs ROC analysis for tripartite rules. If \eqn{plot=TRUE}, the ROC curve is returned.
#' @param Obj An object of class TGST. 
#' @param plot Logical parameter indicating if ROC curve should be plotted. Default is \code{plot=TRUE}. If false, then only AUC is calculated.
#' @return 
#' AUC (the area under ROC curve) and ROC curve.
#' @keywords ROC  
#' @import ggplot2
#' @export
#' @examples
#' d = Simdata
#' Z = d$Z # True Disease Status
#' S = d$S # Risk Score
#' phi = 0.1 #10% of patients taking viral load test
#' lambda = 0.5
#' Obj = TGST(Z, S, phi, method="nonpar")
#' ROCAnalysis(Obj, plot=TRUE)


ROCAnalysis <- function(Obj,plot=TRUE){

  rules <- Obj@Rules
  
  fnr.fpr <- Obj@FNR.FPR
  
  auc <- cal.AUC(Obj@Z,Obj@S,rules[,1],rules[,2])
  
  if(plot==TRUE){  
    #ROC curve
    dat=data.frame(fnr.fpr[,2], 1-fnr.fpr[,1])
    names(dat) <- c("FPR","TPR")

    x = min(dat$FPR)+3/4*(max(dat$FPR)-min(dat$FPR))
    y = min(dat$TPR)+1/4*(max(dat$TPR)-min(dat$TPR))

    dat.add = data.frame(rbind(c(0,0), dat, c(1,1)))
    names(dat.add) = c("FPR","TPR")
    d <- ggplot() +
         geom_line(aes(dat.add$FPR,dat.add$TPR),color="lightblue4", lwd=1) + 
         ylim(0,1)
    
    #d <- ggplot(dat, aes(FPR,TPR)) +geom_ribbon(aes(ymin=0, ymax=TPR), fill="lightblue3", color="lightpink3")+geom_line(color="lightblue4", lwd=1)
    d <- d+geom_text(aes(x,y, label = paste("AUC=",round(auc,3),sep=" ")))    
    plot(d)
  }
  Output = auc
  names(Output) = "AUC"
  return(Output) 
}

