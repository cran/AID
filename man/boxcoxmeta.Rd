\name{boxcoxmeta}
\alias{boxcoxmeta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ensemble Based Box-Cox Transformation via Meta Analysis for Normality of a Variable
}
\description{
\code{boxcoxmeta} performs ensemble based Box-Cox transformation via meta analysis for normality of a variable and provides graphical analysis.  
}
\usage{boxcoxmeta(data, lambda = seq(-3,3,0.01), nboot = 100, lambda2 = NULL, plot = TRUE, 
  alpha = 0.05, verbose = TRUE)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a numeric vector of data values.}
  \item{lambda}{a vector which includes the sequence of candidate lambda values. Default is set to (-3,3) with increment 0.01.} 
  \item{nboot}{a number of Bootstrap samples to estimate standard errors of lambda estimates.}   
  \item{lambda2}{a numeric for an additional shifting parameter. Default is set to lambda2 = 0.} 
  \item{plot}{a logical to plot histogram with its density line and qqplot of raw and transformed data. Defaults plot = TRUE.} 
  \item{alpha}{the level of significance to check the normality after transformation. Default is set to alpha = 0.05.}
  \item{verbose}{a logical for printing output to R console.}

}

\details{

Denote \eqn{y} the variable at the original scale and \eqn{y'} the transformed variable. The Box-Cox power transformation is defined by:
  
  \deqn{y' = \left\{ \begin{array}{ll}
    \frac{y^\lambda - 1}{\lambda} \mbox{ ,  if $\lambda \neq 0$} \cr
    log(y) \mbox{ , if $\lambda = 0$} 
    \end{array} \right.}


  If the data include any nonpositive observations, a shifting parameter \eqn{\lambda_2}{lambda2} can be included in the transformation given by: 
  
  \deqn{y' = \left\{ \begin{array}{ll}
    \frac{(y + \lambda_2)^\lambda - 1}{\lambda} \mbox{ ,  if $\lambda \neq 0$} \cr
    log(y + \lambda_2) \mbox{ , if $\lambda = 0$} 
    \end{array} \right.}



}
\value{
A list with class "boxcoxmeta" containing the following elements:

\item{method}{name of method}
\item{lambda.hat}{estimate of Box-Cox Power transformation parameter}
\item{lambda2}{additional shifting parameter}
\item{result}{a data frame containing the result}
\item{alpha}{the level of significance to assess normality.} 
\item{tf.data}{transformed data set}
\item{var.name}{variable name}


}
\references{
Yilmaz, M. A., Dag, O. (2022). Ensemble Based Box-Cox Transformation via Meta Analysis. \emph{Journal of Advanced Research in Natural and Applied Sciences}, \bold{8:3}, 463--471.
}
\author{
Muhammed Ali Yilmaz, Osman Dag 
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%\code{\link{boxcox}}
%\code{\link{powerTransform}}
%}


\examples{
library(AID)
data(textile)
\donttest{
out <- boxcoxmeta(textile[,1])
out$lambda.hat # the estimate of Box-Cox parameter 
out$tf.data # transformed data set
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{functions}

