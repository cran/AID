\name{AADT}
\alias{AADT}
\docType{data}
\title{
Average Annual Daily Traffic Data
}
\description{
Average annual daily traffic data collected from the Minnesota Department of Transportation data base. 
}
\usage{data(AADT)}
\format{
  A data frame with 121 observations on the following 8 variables.
  \describe{
    \item{\code{aadt}}{average annual daily traffic for a section of road}
    \item{\code{ctypop}}{population of county}
    \item{\code{lanes}}{number of lanes in the section of road}
    \item{\code{width}}{width of the section of road (in feet)}
    \item{\code{control}}{a factor with levels: 1 = access control; 2 = no access control}
    \item{\code{class}}{a factor with levels: 1 = rural interstate; 2 = rural noninterstate; 3 = urban interstate; 4 = urban noninterstate}
    \item{\code{truck}}{availability situation of road section to trucks}
    \item{\code{locale}}{a factor with levels: 1 = rural; 2 = urban, population <= 50,000; 3 = urban, population > 50,000}
  }
}
%%  \details{
%%  ~~ If necessary, more details than the __description__ above %%  ~~
%%  }


\references{

Cheng, C. (1992). Optimal Sampling for Traffic Volume Estimation, Unpublished Ph.D. dissertation, University of Minnesota, Carlson School of Management.

Neter, J., Kutner, M.H., Nachtsheim, C.J.,Wasserman, W. (1996). \emph{Applied Linear Statistical Models (4th ed.)}, Irwin, page 483.

}


\examples{

library(AID)

data(AADT)
attach(AADT)
hist(aadt)
out <- boxcoxfr(aadt, class)
confInt(out)
}
\keyword{datasets}
