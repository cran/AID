\name{boxcoxfr}
\alias{boxcoxfr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Box-Cox Transformation for One-Way Independent Groups Designs
}
\description{
\code{boxcoxfr} performs Box-Cox transformation for one-way independent groups designs. It is useful to use if the normality or/and the homogenity of variance is/are not satisfied while comparing two or more groups. 

    
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

\details{
Maximum likelihood estimation in feasible region (MLEFR) is used while estimating transformation parameter. MLEFR maximizes the likehood function in feasible region constructed by Shapiro-Wilk test and Bartlett's test. After transformation, normality of the data in each group and homogeneity of variance are assessed by Shapiro-Wilk test and Bartlett's test, respectively. 
}


\value{
A list with class "boxcoxfr" containing the following elements:

\item{method}{method applied in the algorithm}
\item{lambda.hat}{the estimated lambda}
\item{shapiro}{a data frame which gives the test results for the normality of groups via Shapiro-Wilk test}
\item{bartlett}{a matrix which returns the test result for the homogenity of variance via Bartlett's test}
\item{alpha}{the level of significance to assess the assumptions.} 
\item{tf.data}{transformed data set}
\item{x}{a factor object which gives the group for the corresponding elements of y}
\item{y.name}{variable name of y}
\item{x.name}{variable name of x}

}

\references{

Dag, O., Ilk, O. An Algorithm for Estimating Box-Cox Transformation Parameter in ANOVA. \emph{Communications in Statistics - Simulation and Computation,  Accepted (June 16, 2016)}.

}


\author{
Osman Dag, Ozlem Ilk
}


\examples{

library(AID)

data(AADT)
attach(AADT)
out <- boxcoxfr(aadt, class)
out$shapiro
out$bartlett
out$tf.data
confInt(out, level = 0.95)


data <- rnorm(120, 10, 1)
factor <- rep(c("X", "Y", "Z"), each = 40)
out <- boxcoxfr(data, factor, lam = seq(-5, 5, 0.01), tau = 0.01, alpha = 0.01)
confInt(out, level = 0.95)


}

\keyword{functions}
