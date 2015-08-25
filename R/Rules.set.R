#' Nonparametric Rules Set
#'
#' This function gives you all possible cutoffs \eqn{[l,u]} for tripartite rules, by applying nonparametric search to the given data. 
#' \deqn{P(S in [l,u]) \le \phi}
#' @param Z True disease status (No disease / treatment success coded as Z=0, diseased / treatment failure coded as Z=1). 
#' @param S Risk score. 
#' @param phi Percentage of patients taking viral load test. 
#' @return 
#' Matrix with 4 columns. Each row is a possible tripartite rule, with output on lower cutoff, upper cutoff and corresponding misclassification rates (FNR, FPR).
#' @keywords Nonparametric, tripartite rules, FNR, FPR.
#' @export
#' @examples
#' d = Simdata
#' Z = d$Z # True Disease Status
#' S = d$S # Risk Score
#' phi = 0.1 #10% of patients taking viral load test
#' Rules.set( Z, S, phi)

nonpar.fnr.fpr <- function(Z,S,l,u){
  if(length(l)!=length(u))   #l is lower cutoff, u is upper cutoff
    cat("***** Warning: Wrong rules set. \n")
  #l is lower cutoff, u is upper cutoff
  n.bounds <-  length(l) #number of all possible rules
  mean.Z <- mean(Z,na.rm=TRUE)
  fnr.fpr <- NULL
  for(i in 1:n.bounds){
    fnr.fpr <- rbind(fnr.fpr, c(mean((S<l[i])*Z,na.rm=TRUE)/mean(Z,na.rm=TRUE),
                                mean((S>u[i])*(1-Z),na.rm=TRUE)/mean(1-Z,na.rm=TRUE)))
  }
  return(fnr.fpr)
}


nonpar.rules <- function(Z,S,phi){
  if(length(Z)!=length(S))
    cat("***** Warning: Disease status and risk score vector lengths do not match. \n")
  data <- cbind(Z,S)
  Z <- data[complete.cases(data),1]
  S <- data[complete.cases(data),2]
  
  Z <- 1*Z #make logical values into {0,1}
  if(phi>1 || phi<0)
    cat("***** Warning: Invalid phi. \n")
  if(cor(Z,S)<0)
    cat("***** Warning: Disease status is negatively associated with risk score.
        Suggest using (-) value for risk score. \n")
  
  # "total.sam.sz" stands for total sample size
  total.sam.sz <- length(Z)
  # total unique sort risk.score
  S.srt <- sort(unique(S))
  # length(S)=645, length(risk.sc.srt)=457
  n.unique.S <- length(S.srt)
  # empirical cdf of S, cum.F=c(0,....,1)
  cum.F <- 0
  for(i in 1:n.unique.S){
    cum.F <- c(cum.F, mean(S<=S.srt[i],na.rm=TRUE))
  }
  # this can be also achieved by 
  # cdf.S=ecdf(S)
  # cum.F = c(0,cdf.S(S.srt))
  
  # identify all cut-off points that allow the percent in the middle
  # to be no more than phi
  
  # by cum.F[i+1]=P(S<=S.srt[i])=G(S.srt[i]) and cum.F[1]=0
  # the following code returns bounds=c(i,j), where
  # cum.F[j+1]-cum.F[i]=G(S.srt[j])-G(S.srt[i-1])
  # =P(S \in ( S.srt[i-1] , S.srt[j] ])=P(S \in [ S.srt[i] , S.srt[j] ])
  # <=\phi 
  flg <- T; i <- 1;j1 <- 0; bounds <- NULL
  while(flg){
    j <- sum((cum.F-cum.F[i]) <=phi)-1
    if(j>=i){
      if(i == 1){
        bounds <- c(i, j) 
        j1 <- j
      }
      if(i>1){
        if(j>j1){
          bounds <- rbind(bounds, c(i,j))
        }
        j1 <- j
      }
    }
    #c(i,j)
    i <- i+1
    if(i>n.unique.S || j1==n.unique.S) flg=F
  }
  if(phi==0){
    bounds <- cbind(1:n.unique.S, 1:n.unique.S)
    cat("***** Warning: 0 patient taking viral load test. \n")
  }
  l <- S.srt[bounds[,1]]
  u <- S.srt[bounds[,2]]
  return(cbind(l,u))
}


Rules.set <- function(Z,S,phi){
  
  rules <- nonpar.rules(Z,S,phi)
  fnr.fpr <- nonpar.fnr.fpr(Z,S,rules[,1],rules[,2])
  outpt <- cbind(rules, fnr.fpr)

  ###
  outpt <- cbind(outpt)
  colnames(outpt) <- c("lower.cutoff", "upper.cutoff", "FNR", "FPR")
  return(outpt)
}

