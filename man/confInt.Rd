\name{confInt}
\alias{confInt}
\alias{confInt.default}
\alias{confInt.boxcoxnc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mean and Confidence Interval for Back Transformed Data
}
\description{
\code{confInt} is a generic function to calculate mean and confidence interval for back transformed data.
}

\usage{
\method{confInt}{boxcoxnc}(x, level = 0.95, verbose = TRUE, ...)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a \code{boxcoxnc} object.}
  \item{level}{the confidence level.} 
  \item{verbose}{a logical for printing output to R console.} 
  \item{...}{additional argument(s) for methods.}

}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}


\value{
Returns a matrix of output.
}



\author{
Osman Dag
}



\examples{

library(AID)

data(textile)
out <- boxcoxnc(textile[,1])
confInt(out) # mean and confidence interval for back transformed data

}

\keyword{functions}
