#' Exploratory Analysis
#'
#' This function allows you to compute the percentage of treatment failure (ART=1) and show distribution summary of risk score by ART status.
#' @param ART Treatment status. 
#' @param S Risk score. 
#' @param Phi Percentage of taking viral load test. 
#' @param RiskF Risk function.
#' @return Percentage of treatment failure (ART=1); histogram and summary of risk score by ART status.
#' @keywords Exploratory analysis, risk score distribution.
#' @export
#' @examples
#' Explore(ART,S=NULL,phi=NULL,RiskF=NULL)

Explore <- function(ART,S=NULL,phi=NULL,RiskF=NULL){
  percF = mean(ART,na.rm=TRUE)
  return(percF)
}