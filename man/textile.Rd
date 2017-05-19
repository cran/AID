\name{textile}
\alias{textile}
\docType{data}
\title{
Textile Data
}
\description{
Number of Cycles to Failure of Worsted Yarn
}
\usage{data(textile)}
\format{
  A data frame with 27 observations on the following variable.
  \describe{
    \item{\code{textile}}{a numeric vector for the number of cycles}
  }
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
\references{
Box, G. E. P., Cox, D. R. (1964). An Analysis of Transformations (with discussion). \emph{Journal of the Royal Statistical Society, Series B (Methodological)}, \bold{26}, 211--252.
}
\examples{

library(AID)

data(textile)
hist(textile[,1])
out <- boxcoxnc(textile[,1])
confInt(out)
}
\keyword{datasets}
