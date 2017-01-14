\name{boxcoxfr}
\alias{boxcoxfr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Function to Estimate Box-Cox Transformation Parameter for Comparing Two or More Groups
}
\description{
\code{boxcoxfr} estimates Box-Cox power transformation parameter. It is useful to use if the normality or/and the homogenity of variance is/are not satisfied while comparing two or more groups. While estimating transformation parameter, maximum likelihood estimation in feasible region (MLEFR) is used. MLEFR maximizes the likehood function in feasible region constructed by Shapiro-Wilk test and Bartlett's test. After transformation parameter is estimated, normality of groups and homogenity of variance are checked by Shapiro-Wilk test and Bartlett's test, respectively. 

Dag, O., Ilk, O. An Algorithm for Estimating Box-Cox Transformation Parameter in ANOVA. \emph{Communications in Statistics - Simulation and Computation,  Accepted (June 16, 2016)}.
    

}

\usage{boxcoxfr(y, x, option = "both", lam = seq(-3, 3, 0.01), tau = 0.05, alpha = 0.05, 
  verbose = TRUE)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{a numeric vector of data values.}
  \item{x}{a vector or factor object which gives the group for the corresponding elements of y.}
  \item{option}{a character string to select the desired option for the objective of transformation. "norm" and "var" are the options which search for a transformation to satisfy the normality of groups and the homogenity of variances, respectively. "both" is the option which searches for a transformation to satisfy both the normality of groups and the homogenity of variances. Default is set to "both".}
  \item{lam}{a vector which includes the sequence of feasible lambda values. Default is set to (-3, 3) with increment 0.01.} 
  \item{tau}{the feasible region parameter for the construction of feasible region. Default is set to 0.05. If tau = 0, it returns the MLE of transformation parameter.} 
  \item{alpha}{the level of significance to check the normality and variance homogenity after transformation. Default is set to alpha = 0.05.}
  \item{verbose}{a logical for printing output to R console.}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}


\value{
Returns a list containing following elements:

\item{method}{method applied in the algorithm}
\item{lambda.hat}{the estimated lambda}
\item{shapiro}{a data frame which gives the test results for the normality of groups via Shapiro-Wilk test}
\item{bartlett}{a matrix which returns the test result for the homogenity of variance via Bartlett's test}
\item{tf.data}{transformed data set}
\item{x}{a factor object which gives the group for the corresponding elements of y}
\item{date}{date at which the function is used}

}

\references{

Bartlett, M.S.(1937). Properties of Sufficiency and Statistical Tests. \emph{Proceedings of the Royal Society of London, Series A}, \bold{160}, 268-282. 

Box, G.E.P., Cox, D.R. (1964). An Analysis of Transformations (with discussion).  \emph{Journal of The Royal Statistical Society, Series B (Methodological)}, \bold{26}, 211-252.

Dag, O., Ilk, O. An Algorithm for Estimating Box-Cox Transformation Parameter in ANOVA. \emph{Communications in Statistics - Simulation and Computation,  Accepted (June 16, 2016)}.

Shapiro, S.S., Wilk, M.B. (1965). An Analysis of Variance Test for Normality (Complete Samples). \emph{Biometrika}, \bold{52(3/4)}, 591-611.

}


\author{
Osman Dag, Ozlem Ilk
}


\examples{

data(AADT)
attach(AADT)
out=boxcoxfr(aadt, class)
out$shapiro
out$bartlett
out$tf.data # transformed data set

data = rnorm(120, 10, 1)
factor = rep(c("X", "Y", "Z"), each = 40)
boxcoxfr(data, factor, lam = seq(-5, 5, 0.01), tau = 0.01, alpha = 0.01)

}

\keyword{functions}
