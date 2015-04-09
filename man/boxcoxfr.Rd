\name{boxcoxfr}
\alias{boxcoxfr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Function to Estimate Box-Cox Transformation Parameter for Comparing Two or More Groups
}
\description{
\code{boxcoxfr} estimates Box-Cox power transformation parameter. It is useful to use if the normality or/and the homogenity of variance is/are not satisfied while comparing two or more groups. While estimating transformation parameter, maximum likelihood estimation in feasible region (MLEFR) is used. MLEFR maximizes the likehood function in feasible region constructed by Shapiro-Wilk test and Bartlett's test. After transformation parameter is estimated, normality of groups and homogenity of variance are checked by Shapiro-Wilk test and Bartlett's test, respectively.     

}

\usage{boxcoxfr(y, x, option = "both", lam = seq(-2, 2, 0.02), alpha = 0.05)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{a numeric vector of data values.}
  \item{x}{a vector or factor object which gives the group for the corresponding elements of y.}
  \item{option}{a character string to select the desired option for the objective of transformation. "norm" and "var" are the options which search for a transformation to satisfy the normality of groups and the homogenity of variances, respectively. "both" is the option which searches for a transformation to satisfy both the normality of groups and the homogenity of variances. Default is set to "both".}
  \item{lam}{a vector which includes the sequence of feasible lambda values. Default is set to (-2, 2) with increment 0.02.} 
  \item{alpha}{the significance level for the construction of feasible region.} 

}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}


\value{
Returns a list containing following elements:
\item{method}{method applied in the algorithm}
\item{date}{date which the function is used}
\item{lambda.hat}{the estimated lambda}
\item{shapiro.test}{a matrix which gives the test results for the normality of groups via shapiro.test}
\item{bartlett.test}{a matrix which returns the test result for the homogenity of variance via bartlett.test}
}

\references{

Bartlett, M.S.(1937). Properties of Sufficiency and Statistical Tests. \emph{Proceedings of the Royal Society of London, Series A}, \bold{160}, 268-282. 

Box, G.E.P., Cox, D.R. (1964). An Analysis of Transformations (with discussion).  \emph{Journal of The Royal Statistical Society, Series B (Methodological)}, \bold{26}, 211-252.

Shapiro, S.S., Wilk, M.B. (1965). An Analysis of Variance Test for Normality (Complete Samples). \emph{Biometrika}, \bold{52(3/4)}, 591-611.

}


\author{
Osman Dag, Ozlem Ilk
}


\examples{

data = rnorm(120, 10, 1)
factor = rep(c("X", "Y", "Z"), each = 40)
boxcoxfr(data, factor, lam = seq(-3, 3, 0.01))

}

\keyword{functions}
