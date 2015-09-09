#' Create a TVLT Object
#'
#' Create a TVLT object, usually used as an input for optimal rule search and ROC analysis. 
#' @param Z A vector of true disease status (Viral failure coded as Z=1).
#' @param S A vector of risk Score.
#' @param phi Percentage of patients taking viral load test. 
#' @param method Method for searching for the optimal tripartite rule, options are "nonpar" (default) and "semipar".
#' @return 
#' An object of class \code{TVLT}.The class contains 6 slots: phi (percentage of VL tests), Z (true viral failure status), S (risk score), Rules (all possible tripartite rules), Nonparametric (logical indicator of the approach), and FNR.FPR (miclassification rates).
#' @keywords Tripartite rules, FNR, FPR.
#' @export
#' @examples
#' d = Simdata
#' Z = d$Z # True Disease Status
#' S = d$S # Risk Score
#' phi = 0.1 #10% of patients taking viral load test
#' TVLT( Z, S, phi, method="nonpar")
#' 

TVLT <- function(Z, S, phi, method="nonpar"){
  
  #consider only complete cases
  data <- cbind(Z,S)
  Z <- data[complete.cases(data),1]
  S <- data[complete.cases(data),2]
  
  if(!(method %in% c("nonpar","semipar"))) stop("***** Warning: Wrong method input. \n")

  Z <- Z*1
  rules <- nonpar.rules(Z,S,phi)
  if(method=="nonpar"){
    logic <- TRUE
    fnr.fpr <- nonpar.fnr.fpr(Z,S,rules[,1],rules[,2])
  }
  if(method=="semipar"){
    logic <- FALSE
    fnr.fpr <-  semipar.fnr.fpr(Z,S,rules[,1],rules[,2])
  }
  colnames(fnr.fpr) <- c("FNR","FPR")
  result=new("TVLT", phi=phi, Z=Z, S=S, Rules=rules, Nonparametric=logic, FNR.FPR=fnr.fpr)
  
  return(result)
}

