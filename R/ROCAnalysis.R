#' ROC Analysis
#'
#' This function performs ROC analysis for tripartite rules. If \eqn{plot=TRUE}, the ROC curve is returned.
#' @param Obj An object of class TGST. 
#' @param plot Logical parameter indicating if ROC curve should be plotted. Default is \code{plot=TRUE}. If false, then only AUC is calculated.
#' @return 
#' AUC (the area under ROC curve) and ROC curve.
#' @keywords ROC, AUC.
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
    plot(fnr.fpr[,2],1-fnr.fpr[,1],type="l",xlab="FPR",ylab="TPR",main="ROC Curve")
    legend('bottomright',paste("AUC=",round(auc,3),sep=" "),bty ="n",cex=0.8)
  }
  Output = auc
  names(Output) = "AUC"
  invisible(Output)  
}

