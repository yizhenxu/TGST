## ----eval=FALSE---------------------------------------------------------------
#  d = Simdata
#  Z = d$Z # True Disease Status
#  S = d$S # Risk Score
#  Check.exp.tilt( Z, S)

## ----eval=FALSE---------------------------------------------------------------
#  d = Simdata
#  Z = d$Z # True Disease Status
#  S = d$S # Risk Score
#  phi = 0.1 #10\% of patients taking viral load test
#  lambda =0.5
#  CV.TVLT(Z, S, phi, K = 10, method = "semipar", lambda)

## ----eval=FALSE---------------------------------------------------------------
#  d = Simdata
#  Z = d$Z # True Disease Status
#  S = d$S # Risk Score
#  phi = 0.1 #10\% of patients taking viral load test
#  lambda = 0.5
#  Opt.nonpar.rule( Z, S, phi, lambda)

## ----eval=FALSE---------------------------------------------------------------
#  d = Simdata
#  Z = d$Z # True Disease Status
#  S = d$S # Risk Score
#  phi = 0.1 #10% of patients taking viral load test
#  lambda = 0.5
#  Opt.semipar.rule( Z, S, phi, lambda)

## ----eval=FALSE---------------------------------------------------------------
#  d = Simdata
#  Z = d$Z # True Disease Status
#  S = d$S # Risk Score
#  phi = 0.1 #10% of patients taking viral load test
#  a = ROC.nonpar( Z, S, phi,plot=TRUE)
#  a$AUC
#  a$FNR
#  a$FPR

## ----eval=FALSE---------------------------------------------------------------
#  d = Simdata
#  Z = d$Z # True Disease Status
#  S = d$S # Risk Score
#  phi = 0.1 #10% of patients taking viral load test
#  a = ROC.semipar( Z, S, phi,plot=TRUE)
#  a$AUC
#  a$FNR
#  a$FPR

## ----eval=FALSE---------------------------------------------------------------
#  d = Simdata
#  Z = d$Z # True Disease Status
#  S = d$S # Risk Score
#  phi = 0.1 #10\% of patients taking viral load test
#  Rules.set( Z, S, phi)

## ----eval=FALSE---------------------------------------------------------------
#  d = Simdata
#  Z = d$Z # True Disease Status
#  S = d$S # Risk Score
#  phi = 0.1 #10\% of patients taking viral load test
#  Semi.par.rule( Z, S, phi)

## ----eval=FALSE---------------------------------------------------------------
#  data(Simdata)
#  summary(Simdata)
#  plot(Simdata)

## ----eval=FALSE---------------------------------------------------------------
#  d = Simdata
#  Z = d$Z # True Disease Status
#  S = d$S # Risk Score
#  TVLT.summ(Z,S)

