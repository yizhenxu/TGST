#' Nonparametric Rules
#'
#' This function gives you all possible cutoffs \eqn{[l,u]} for tripartite rules, by applying nonparametric search to the given data. 
#' \deqn{P(S in [l,u]) \le \phi}
#' @param Z True disease status (No disease / treatment success coded as Z=0, diseased / treatment failure coded as Z=1). 
#' @param S Risk score. 
#' @param phi Percentage of patients taking viral load test. 
#' @return 
#' Cutoff points l and u;
#' Misdiagnoses rate for viral failure (i.e., false negative rate, FNR) and otherwise (i.e., false positive rate, FPR). 
#' @keywords Nonparametric, tripartite rules, FNR, FPR.
#' @export
#' @examples
#' data = Simdata
#' Z = d$Z # True Disease Status
#' S = d$S # Risk Score
#' phi = 0.1 #10% of patients taking viral load test
#' Rules.set( Z, S, phi)

Rules.set <- function(Z,S,phi){
  Z <- 1*Z #make logical values into {0,1}
  if(cor(Z,S,use = "complete.obs")<0)
    cat("***** Warning: Disease status is negatively associated with risk score.
        Suggest using (-) value for risk score. \n")
  # "total.sam.sz" stands for total sample size
  total.sam.sz <- length(Z)
  # "vl.test.sam.sz" stands for the segments size of the sample
  vl.test.sam.sz <- floor(total.sam.sz*phi)
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
  flg <- T; i <- 1
  while(flg){
    j <- sum((cum.F-cum.F[i]) <=phi)-1
    if(i == 1){
      bounds <- c(i, j) #j=63
      j1 <- j
    }
    if(i>1){
      if(j>j1){
        bounds <- rbind(bounds, c(i,j))
      }
      j1 <- j
    }
    #c(i,j)
    i <- i+1
    if(i>n.unique.S || j1==n.unique.S) flg=F
  }
  if(phi==0){
    bounds <- cbind(1:n.unique.S, 1:n.unique.S)
    cat("***** Warning: 0 patient taking viral load test. \n")
  }
  # FNR and FPR
  n.bounds <- dim(bounds)[1] #number of bounds
  mean.Z <- mean(Z,na.rm=TRUE)
  fnr.fpr <- NULL
  for(i in 1:n.bounds){
    fnr.fpr <- rbind(fnr.fpr, c(mean((S<S.srt[bounds[i,1]])*Z,na.rm=TRUE)/mean(Z,na.rm=TRUE),
                                mean((S>S.srt[bounds[i,2]])*(1-Z),na.rm=TRUE)/mean(1-Z,na.rm=TRUE)))
  }
  outpt <- cbind(S.srt[bounds[,1]], S.srt[bounds[,2]], fnr.fpr)

  ###
  outpt <- cbind(outpt)
  colnames(outpt) <- c("lower.cutoff", "upper.cutoff", "FNR", "FPR")
  return(outpt)
}

