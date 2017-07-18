\name{boxcoxlm}
\alias{boxcoxlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Box-Cox Transformation for Linear Models
}
\description{
\code{boxcoxlm} performs Box-Cox transformation for linear models and provides graphical analysis of residuals after transformation.  
}
\usage{boxcoxlm(x, y, method = "lse", lambda = seq(-3,3,0.01), lambda2 = NULL, plot = TRUE, 
  alpha = 0.05, verbose = TRUE)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a nxp matrix, n is the number of observations and p is the number of variables.}
  \item{y}{a vector of response variable.}
  \item{method}{a character string to select the desired method to be used to estimate Box-Cox transformation parameter. To use Shapiro-Wilk test method should be set to "sw". For method = "ad", boxcoxnc function uses Anderson-Darling test to estimate Box-Cox transformation parameter. Similarly, method should be set to "cvm", "pt", "sf", "lt", "jb", "mle", "lse" to use Cramer-von Mises, Pearson Chi-square, Shapiro-Francia, Lilliefors and Jarque-Bera tests, maximum likelihood estimation and least square estimation, respectively. Default is set to method = "lse".}
  \item{lambda}{a vector which includes the sequence of candidate lambda values. Default is set to (-3,3) with increment 0.01.} 
  \item{lambda2}{a numeric for an additional shifting parameter. Default is set to lambda2 = 0.} 
  \item{plot}{a logical to plot histogram with its density line and qqplot of residuals before and after transformation. Defaults plot = TRUE.} 
  \item{alpha}{the level of significance to assess the normality of residuals after transformation. Default is set to alpha = 0.05.}
  \item{verbose}{a logical for printing output to R console.}

}

\details{

Denote \eqn{y} the variable at the original scale and \eqn{y'} the transformed variable. The Box-Cox power transformation is defined by:
  
  \deqn{y' = \left\{ \begin{array}{ll}
    \frac{y^\lambda - 1}{\lambda} = \beta_0 + \beta_1x_1 + ... + \epsilon \mbox{ ,  if $\lambda \neq 0$} \cr
    log(y) = \beta_0 + \beta_1x_1 + ... + \epsilon \mbox{ , if $\lambda = 0$} 
    \end{array} \right.}


  If the data include any nonpositive observations, a shifting parameter \eqn{\lambda_2}{lambda2} can be included in the transformation given by: 
  
  \deqn{y' = \left\{ \begin{array}{ll}
    \frac{(y + \lambda_2)^\lambda - 1}{\lambda} = \beta_0 + \beta_1x_1 + ... + \epsilon \mbox{ ,  if $\lambda \neq 0$} \cr
    log(y + \lambda_2) = \beta_0 + \beta_1x_1 + ... + \epsilon \mbox{ , if $\lambda = 0$} 
    \end{array} \right.}

Maximum likelihood estimation and least square estimation are equivalent while estimating Box-Cox power transformation parameter (Kutner et al., 2005). Therefore, these two methods return the same result.

}
\value{
A list with class "boxcoxlm" containing the following elements:

\item{method}{method preferred to estimate Box-Cox transformation parameter}
\item{lambda.hat}{estimate of Box-Cox Power transformation parameter based on corresponding method}
\item{lambda2}{additional shifting parameter}
\item{statistic}{statistic of normality test for residuals after transformation based on specified normality test in method. For mle and lse, statistic is obtained by Shapiro-Wilk test for residuals after transformation}
\item{p.value}{p.value of normality test for residuals after transformation based on specified normality test in method. For mle and lse, p.value is obtained by Shapiro-Wilk test for residuals after transformation}
\item{alpha}{the level of significance to assess normality of residuals} 
\item{tf.y}{transformed response variable}
\item{tf.residuals}{residuals after transformation}
\item{y.name}{response name}
\item{x.name}{x matrix name}

}
\references{
Asar, O., Ilk, O., Dag, O. (2017). Estimating Box-Cox Power Transformation Parameter via Goodness of Fit Tests. \emph{Communications in Statistics - Simulation and Computation}, \bold{46:1}, 91--105.

Kutner, M. H., Nachtsheim, C., Neter, J., Li, W. (2005). \emph{Applied Linear Statistical Models. (5th ed.)}. New York: McGraw-Hill Irwin.

}
\author{
Osman Dag, Ozlem Ilk
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%\code{\link{boxcox}}
%\code{\link{powerTransform}}
%}
\examples{

library(AID)

trees=as.matrix(trees)
boxcoxlm(x = trees[,1:2], y = trees[,3])

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{functions}

