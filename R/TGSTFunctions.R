
#All possible rules

#' Nonparametric Rules Set
#'
#' This function gives you all possible cutoffs \eqn{[l,u]} for tripartite rules, by applying nonparametric search to the given data. 
#' \deqn{P(S in [l,u]) \le \phi}
#' @param Z True disease status (No disease / treatment success coded as Z=0, diseased / treatment failure coded as Z=1). 
#' @param S Risk score. 
#' @param phi Percentage of patients taking viral load test. 
#' @return 
#' Matrix with 2 columns. Each row is a possible tripartite rule, with output on lower and upper cutoff.
#' @keywords nonparametric rules
#' @export
#' @examples
#' d = Simdata
#' Z = d$Z # True Disease Status
#' S = d$S # Risk Score
#' phi = 0.1 #10% of patients taking viral load test
#' nonpar.rules( Z, S, phi)

nonpar.rules <- function(Z,S,phi){
  if(length(Z)!=length(S))
    message("***** Warning: Disease status and risk score vector lengths do not match. \n")
  data <- cbind(Z,S)
  Z <- data[complete.cases(data),1]
  S <- data[complete.cases(data),2]
  
  Z <- 1*Z #make logical values into {0,1}
  if(phi>1 || phi<0)
    message("***** Warning: Invalid phi. \n")
  if(cor(Z,S)<0)
    message("***** Warning: Disease status is negatively associated with risk score.
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
    message("***** Warning: 0 patient taking viral load test. \n")
  }
  l <- S.srt[bounds[,1]]
  u <- S.srt[bounds[,2]]
  return(cbind(l,u))
}




#' Nonparametric FNR FPR of the rules
#'
#' This function gives you the nonparametric FNR and FPR associated with a given tripartite rule.
#' @param Z True disease status (No disease / treatment success coded as Z=0, diseased / treatment failure coded as Z=1). 
#' @param S Risk score. 
#' @param l Lower cutoff of tripartite rule. 
#' @param u Upper cutoff of tripartite rule. 
#' @return 
#' Matrix with 2 columns. Each row is a set of nonparametric (FNR, FPR) on an associated tripartite rule.
#' @keywords nonparametric rules
#' @export
#' @examples
#' d = Simdata
#' Z = d$Z # True Disease Status
#' S = d$S # Risk Score
#' phi = 0.1 #10% of patients taking viral load test
#' rules = nonpar.rules( Z, S, phi)
#' nonpar.fnr.fpr(Z,S,rules[1,1],rules[1,2])

nonpar.fnr.fpr <- function(Z,S,l,u){
  if(length(l)!=length(u))   #l is lower cutoff, u is upper cutoff
    message("***** Warning: Wrong rules set. \n")
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



#' Semiparametric FNR FPR of the rules
#'
#' This function gives you the semiparametric FNR and FPR associated with a given tripartite rule.
#' @param Z True disease status (No disease / treatment success coded as Z=0, diseased / treatment failure coded as Z=1). 
#' @param S Risk score. 
#' @param l Lower cutoff of tripartite rule. 
#' @param u Upper cutoff of tripartite rule. 
#' @return 
#' Matrix with 2 columns. Each row is a set of semiparametric (FNR, FPR) on an associated tripartite rule.
#' @keywords semiparametric rules
#' @export
#' @examples
#' d = Simdata
#' Z = d$Z # True Disease Status
#' S = d$S # Risk Score
#' phi = 0.1 #10% of patients taking viral load test
#' rules = nonpar.rules( Z, S, phi)
#' semipar.fnr.fpr(Z,S,rules[1,1],rules[1,2])

semipar.fnr.fpr <- function(Z,S,l,u){
  if(length(l)!=length(u))   #l is lower cutoff, u is upper cutoff
    message("***** Warning: Wrong rules set. \n")
  p <- mean(Z)
  temp <- density(S) #marginal density (S)
  fit <- glm(Z~ S, family=binomial)
  beta0star <- fit$coef[1]-log(p/(1-p))
  t <- exp(beta0star+temp$x*fit$coef[2]) #g1=t*g0 under exp tilt assumption
  g1 <- temp$y/(p+(1-p)/t)
  g0 <- temp$y/(p*t+1-p)
  x <- temp$x
  
  len <- length(x)
  dif <- x[2:len]-x[1:(len-1)]
  
  cal.fnr <- function(dens,a){  
    if( a>max(x) ){
      area <- 1
    } else if( a<min(x) ){
      area <- 0
    } else {
      diff <- a-x
      diff1 <- diff[diff<=0][1] 
      indx <- which(diff==diff1)#return index of nearest right endpoint
      area <- sum(dens[1:(indx-2)]*dif[1:(indx-2)])+dens[indx-1]*(a-x[indx-1])
    }
    return(area)
  } 
  
  fnr.fpr <- NULL
  K <- length(l)
  for( i in 1:K){
    fnr <- cal.fnr(g1,l[i])
    fpr <- 1-cal.fnr(g0,u[i])
    fnr.fpr <- rbind(fnr.fpr,c(fnr,fpr))        
  }
  
  return(fnr.fpr)
  
}




#' Calculate AUC
#'
#' This function gives you the AUC associated with the rules set.
#' @param Z True disease status (No disease / treatment success coded as Z=0, diseased / treatment failure coded as Z=1). 
#' @param S Risk score. 
#' @param l Lower cutoff of all possible tripartite rules. 
#' @param u Upper cutoff of all possible tripartite rules. 
#' @return 
#' AUC.
#' @keywords AUC
#' @export
#' @examples
#' d = Simdata
#' Z = d$Z # True Disease Status
#' S = d$S # Risk Score
#' phi = 0.1 #10% of patients taking viral load test
#' rules = nonpar.rules( Z, S, phi)
#' cal.AUC(Z,S,rules[,1],rules[,2])

cal.AUC <- function(Z,S,l,u){
  ## AUC
  #Write the kth rule in Rule.set as (i_k,j_k), let j_0=0
  #Hphi(u)=argmin_w {G(u)-G(w)<=phi}
  #For Sj in (j_{k-1},j_k], Hphi(Sj)=i_k
  n = length(Z)
  p <- mean(Z)
  Hphi <- function(Sj,bounds=cbind(l,u)){
    diff <- Sj-bounds[,2]
    diff1 <- diff[diff<=0][1] #Sj-j_k, where Sj in (j_{k-1},j_k]
    indx <- which(diff==diff1)
    return(bounds[indx,1])
  }
  #calculate AUC from eqn (10) pg1177
  auc <- 0
  for(j in 1:n){
    auc <- auc+sum(Z*(1-Z[j])*((S>Hphi(S[j]))+(S==Hphi(S[j]))/2))
  }
  auc <- auc/(n^2*p*(1-p))  
}