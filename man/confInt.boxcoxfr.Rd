\name{confInt.boxcoxfr}
\alias{confInt.boxcoxfr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mean and Asymmetric Confidence Interval for Back Transformed Data
}
\description{
\code{confInt.boxcoxfr} calculates mean and asymmetric confidence interval for back transformed data in each group and plots their error bars with confidence intervals.
}

\usage{
\method{confInt}{boxcoxfr}(x, level = 0.95, plot = TRUE, xlab = NULL, ylab = NULL, title = NULL, 
  width = NULL, verbose = TRUE, ...)}


%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a \code{boxcoxfr} object.}
  \item{level}{the confidence level.} 
  \item{plot}{a logical to plot error bars with confidence intervals.} 
  \item{xlab}{a label for the x axis, defaults to a description of x.} 
  \item{ylab}{a label for the y axis, defaults to a description of y.} 
  \item{title}{a main title for the plot.} 
  \item{width}{a numeric giving the width of the little lines at the tops and bottoms of the error bars (defaults to 0.15).} 
  \item{verbose}{a logical for printing output to R console.} 
  \item{...}{additional argument(s) for methods.}

}

\details{
Confidence interval in each group is constructed separately.
}


\value{
A matrix with columns giving mean, lower and upper confidence limits for back transformed data. These will be labelled as (1 - level)/2 and 1 - (1 - level)/2 in \% (by default 2.5\% and 97.5\%).
}



\author{
Osman Dag}



\examples{

library(AID)

data(AADT)
attach(AADT)
out <- boxcoxfr(aadt, class)
confInt(out, level = 0.95)
}

\keyword{functions}
