\name{boxcoxnc}
\alias{boxcoxnc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Box-Cox Transformation for Normality of a Variable
}
\description{
\code{boxcoxnc} performs Box-Cox transformation for normality of a variable and provides graphical analysis.  
}
\usage{boxcoxnc(data, method = "sw", lambda = seq(-3,3,0.01), lambda2 = NULL, plot = TRUE, 
  alpha = 0.05, verbose = TRUE)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a numeric vector of data values.}
  \item{method}{a character string to select the desired method to be used to estimate Box-Cox transformation parameter. To use Shapiro-Wilk test method should be set to "sw". For method = "ad", boxcoxnc function uses Anderson-Darling test to estimate Box-Cox transformation parameter. Similarly, method should be set to "cvm", "pt", "sf", "lt", "jb", "ac", "mle" to use Cramer-von Mises, Pearson Chi-square, Shapiro-Francia, Lilliefors, Jarque-Bera tests, artificial covariate method and maximum likelihood estimation, respectively. Default is set to method = "sw".}
  \item{lambda}{a vector which includes the sequence of candidate lambda values. Default is set to (-3,3) with increment 0.01.} 
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
A list with class "boxcoxnc" containing the following elements:

\item{method}{method preferred to estimate Box-Cox transformation parameter}
\item{lambda.hat}{estimate of Box-Cox Power transformation parameter based on corresponding method}
\item{lambda2}{additional shifting parameter}
\item{statistic}{statistic of normality test for transformed data based on specified normality test in method. For artificial covariate method, statistic is obtained by Shapiro-Wilk test for transformed data}
\item{p.value}{p.value of normality test for transformed data based on specified normality test in method. For artificial covariate method, p.value is obtained by Shapiro-Wilk test for transformed data}
\item{alpha}{the level of significance to assess normality.} 
\item{tf.data}{transformed data set}
\item{var.name}{variable name}


}
\references{
Asar, O., Ilk, O., Dag, O. (2017). Estimating Box-Cox Power Transformation Parameter via Goodness of Fit Tests. \emph{Communications in Statistics - Simulation and Computation}, \bold{46:1}, 91--105.

Dag, O., Asar, O., Ilk, O. (2014). A Methodology to Implement Box-Cox Transformation When No Covariate is Available. \emph{Communications in Statistics - Simulation and Computation}, \bold{43:7}, 1740--1759.

}
\author{
Osman Dag, Ozgur Asar, Ozlem Ilk
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%\code{\link{boxcox}}
%\code{\link{powerTransform}}
%}
\examples{

library(AID)

data(textile)

out <- boxcoxnc(textile[,1], method = "sw")
out$lambda.hat # the estimate of Box-Cox parameter based on Shapiro-Wilk test statistic 
out$p.value # p.value of Shapiro-Wilk test for transformed data 
out$tf.data # transformed data set
confInt(out) # mean and confidence interval for back transformed data


out2 <- boxcoxnc(textile[,1], method = "sf")
out2$lambda.hat # the estimate of Box-Cox parameter based on Shapiro-Francia test statistic
out2$p.value # p.value of Shapiro-Francia test for transformed data 
out2$tf.data 
confInt(out2) 
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{functions}

