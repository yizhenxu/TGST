###################################################################################
##' Constructor of TVLT class
##'
##' This class contains results of a run to search for all possible rules. 
##'
##' \describe{
##'   \item{phi}{Percentage of patients taking viral load test.}
##'   \item{Z}{A vector of true disease status (Viral failure coded as Z=1).}
##'   \item{S}{A vector of risk Score.}
##'   \item{Rules}{A matrix of all possible tripartite rules (two cutoffs) derived from the training data set.}
##'   \item{Nonparametric}{A boolean indicating if nonparametric approach should be used in calculating the misclassfication rates. If FALSE, semiparametric approach would be used.}
##'   \item{FNR.FPR}{A matrix with two columns of misclassification rates, FNR and FPR.}
##' }
##'
##' @details
##' If res is the result of rankclust(), each slot of results can be reached by res[k]@@slotname, where
##' k is the number of clusters and slotname is the name of the slot we want to reach (see \link{Output-class}).
##' For the slots ll, bic, icl, res["slotname"] returns a vector of size K containing the values of the
##' slot for each number of clusters.
##'
##'
##' @name TVLT-class
##' @rdname TVLT-class
## @exportClass TVLT
##'
##'
setClass(
  Class="TVLT",
  representation=representation(
    phi="numeric",
    Z="numeric",
    S="numeric",
    Rules="matrix",
    Nonparametric="logical",
    FNR.FPR="matrix"
  ),
  prototype=prototype(
    phi=numeric(0),
    Z=numeric(0),
    S=numeric(0),
    Rules=matrix(nrow=0,ncol=0),
    Nonparametric=TRUE,
    FNR.FPR=matrix(nrow=0,ncol=0)
  )
)


#'
#' summary function.
#' 
#' This function This function gives the summary of the data from \code{TVLT}.
#' 
#' @param object Output object from \code{\link{TVLT}}.
#' @return 
#' Percentage of treatment failure; 
#' Summary statistics (mean, standard deviation, minimum, median, maximum and IQR) of risk score by true disease status; 
#' Distribution plot.
#' @name summary
#' @rdname summary-methods
#' @docType methods
#' @exportMethod summary
#' @aliases summary summary, TVLT-method
setMethod(
  f="summary",
  signature = "TVLT",
  definition = function(object,...) {
    Z = object@Z
    S = object@S
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
    dat=data.frame(Z=as.factor(Z), S=S)
    #library(ggplot2)
    ggplot(dat, aes(x=S, fill=Z)) + geom_density(alpha=.3)+ scale_fill_discrete(name="Viral Failure",
                                                                                breaks=c("0","1"),
                                                                                labels=c("Z=0","Z=1"))
    
    #output
    z = list(Percent_of_Viral_Failure=percF,SummaryS0=summ.S0,SummaryS1=summ.S1)
    return(z)
  }
)



