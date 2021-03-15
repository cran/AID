\name{confInt.boxcoxmeta}
\alias{confInt.boxcoxmeta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mean and Asymmetric Confidence Interval for Back Transformed Data
}
\description{
\code{confInt.boxcoxmeta} calculates mean and asymmetric confidence interval for back transformed data.
}

\usage{
\method{confInt}{boxcoxmeta}(x, level = 0.95, verbose = TRUE, ...)}


%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a \code{boxcoxmeta} object.}
  \item{level}{the confidence level.} 
  \item{verbose}{a logical for printing output to R console.} 
  \item{...}{additional argument(s) for methods.}

}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}



\value{
A matrix with columns giving mean, lower and upper confidence limits for back transformed data. These will be labelled as (1 - level)/2 and 1 - (1 - level)/2 in \% (by default 2.5\% and 97.5\%).
}




\author{
Osman Dag, Muhammed Ali Yilmaz}


\examples{
library(AID)
data(textile)
\donttest{
out <- boxcoxmeta(textile[,1])
confInt(out) # mean and confidence interval for back transformed data
}
}

\keyword{functions}
