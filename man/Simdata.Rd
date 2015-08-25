\name{Simdata}
\alias{Simdata}
\docType{data}
\title{
Simulated data for package illustration
}
\description{
A simulated dataset containing true disease status, CD4 count, and risk score. Simulation parameters are determined by fitting a dataset from  Miriam Immunology Clinic. See details.
}
\usage{data(Simdata)}
\format{
  A data frame with 5000 simulated observations on the following 3 variables.
  \describe{
    \item{\code{Z}}{True disease status, with diseased coded as \eqn{Z=1}.}
    \item{\code{CD4}}{Simulated CD4 count with different distribution with respect to different disease status.}
    \item{\code{S}}{Risk score defined as the negative value of CD4 count.}
  }
}
\details{
We first simulate viral status \eqn{Z} assuming \eqn{Z\sim Bernoulli(p)} with \eqn{p=0.25} which was estimated from the Miriam Immunology Clinic data; and then conditonal on \eqn{Z}, simulate \eqn{{CD4|Z=z}=ceiling(W)} with \eqn{W\sim Gamma(\eta_z,\kappa_z)} where \eqn{\eta} and \eqn{\kappa} are shape and scale parameters.\eqn{(\eta0,\kappa0)=(3.2,152)} and \eqn{(\eta1,\kappa1)=(2.3,133). Parameter values are chosen as maximum likelihood estimates obtained by fitting gamma distributions to the Miriam Immunology Clinic data.} 

}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(Simdata)
## maybe str(Simdata) ; plot(Simdata) ...
}
\keyword{datasets}
