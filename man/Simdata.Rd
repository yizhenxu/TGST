\name{Simdata}
\alias{Simdata}
\docType{data}
\title{
Simulated data for package illustration
}
\description{
A simulated dataset containing true disease status and risk score. See details for simulation setting.
}
\usage{data(Simdata)}
\format{
  A data frame with 8000 simulated observations on the following 2 variables.
  \describe{
    \item{\code{Z}}{True disease status (No disease / treatment success coded as Z=0, diseased / treatment failure coded as Z=1).}
    \item{\code{S}}{Risk score. Higher risk score indicates larger tendency of diseased / treatment failure.}
  }
}
\details{
We first simulate viral status \eqn{Z} assuming \eqn{Z\sim Bernoulli(p)} with \eqn{p=0.25}; and then conditonal on \eqn{Z}, simulate \eqn{{S|Z=z}=ceiling(W)} with \eqn{W\sim Gamma(\eta_z,\kappa_z)} where \eqn{\eta} and \eqn{\kappa} are shape and scale parameters.\eqn{(\eta0,\kappa0)=(2.3,80)} and \eqn{(\eta1,\kappa1)=(9.2,62)}. 
}
\examples{
data(Simdata)
## maybe str(Simdata) ; plot(Simdata) ...
}
\keyword{datasets}
