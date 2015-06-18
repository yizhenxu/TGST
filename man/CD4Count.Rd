\name{CD4Count}
\alias{CD4Count}
\docType{data}
\title{
Viral load, CD4 count and CD4 percent of 645 patients.
}
\description{
A data set containing viral load count, CD4 count, CD4 percent, CD4 count 6 months ago and CD4 percent 6 months ago.
}
\usage{data("CD4Count")}
\format{
  A data frame with 645 observations on the following 6 variables.
  \describe{
    \item{\code{PtPF}}{Patient ID}
    \item{\code{PVL_Count_Now}}{Viral load count}
    \item{\code{CD4_count_Now}}{CD4 count at current visit}
    \item{\code{CD4Percent_Now}}{CD4 percent at current visit}
    \item{\code{CD4Count_6months_ago}}{CD4 count six months ago}
    \item{\code{CD4perc_6months_ago}}{CD4 percent six months ago}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(CD4Count)
## maybe str(CD4Count) ; plot(CD4Count) ...
}
\keyword{datasets}
