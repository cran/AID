\name{boxcoxnc}
\alias{boxcoxnc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Box-Cox Transformation for Normality of a Univariate Variable
}
\description{
\code{boxcoxnc} utilizes seven different normality tests and artificial covariate method to estimate Box-Cox power transformation parameter and provides graphical analysis

Asar, O., Ilk, O., Dag, O. (2017). Estimating Box-Cox Power Transformation Parameter via Goodness of Fit Tests. \emph{Communications in Statistics - Simulation and Computation}, \bold{46:1}, 91--105.

Dag, O., Asar, O., Ilk, O. (2014). A Methodology to Implement Box-Cox Transformation When No Covariate is Available. \emph{Communications in Statistics - Simulation and Computation}, \bold{43:7}, 1740--1759.

}
\usage{boxcoxnc(data, method = "sw", lam = seq(-3,3,0.01), plot = TRUE, alpha = 0.05, 
  verbose = TRUE)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a numeric vector of data values.}
  \item{method}{a character string to select the desired method to be used to estimate Box-Cox transformation parameter. To use Shapiro-Wilk test method should be set to "sw". For method = "ad", boxcoxnc function uses Anderson-Darling test to estimate Box-Cox transformation parameter. Similarly, method should be set to "cvm", "pt", "sf", "lt", "jb", "ac" to use Cramer-von Mises, Pearson Chi-square, Shapiro-Francia, Lilliefors, Jarque-Bera tests and artificial covariate method, respectively. Default is set to method = "sw".}
  \item{lam}{a vector which includes the sequence of candidate lambda values. Default is set to (-3,3) with increment 0.01.} 
  \item{plot}{plots histogram and qqplot of raw and transformed data. Defaults plot = TRUE.} 
  \item{alpha}{the level of significance to check the normality after transformation. Default is set to alpha = 0.05.}
  \item{verbose}{a logical for printing output to R console.}

}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Returns a list containing following elements:

\item{method}{method preferred to estimate Box-Cox transformation parameter}
\item{lambda.hat}{estimate of Box-Cox Power transformation parameter based on corresponding method}
\item{statistic}{statistic of normality test for transformed data based on specified normality test in method. For artificial covariate method, statistic is obtained by Shapiro-Wilk test for transformed data}
\item{p.value}{p.value of normality test for transformed data based on specified normality test in method. For artificial covariate method, p.value is obtained by Shapiro-Wilk test for transformed data}
\item{tf.data}{transformed data set}


}
\references{
Asar, O., Ilk, O., Dag, O. (2017). Estimating Box-Cox Power Transformation Parameter via Goodness of Fit Tests. \emph{Communications in Statistics - Simulation and Computation}, \bold{46:1}, 91--105.

Box, G. E. P., Cox, D. R. (1964). An Analysis of Transformations (with discussion). \emph{Journal of the Royal Statistical Society, Series B (Methodological)}, \bold{26}, 211--252.

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
data(textile)

out <- boxcoxnc(textile[,1], method = "sw")
out$lambda.hat # the estimate of Box-Cox parameter based on Shapiro-Wilk test statistic 
out$p.value # p.value of Shapiro-Wilk test for transformed data 
out$tf.data # transformed data set

out2 <- boxcoxnc(textile[,1], method = "sf")
out2$lambda.hat # the estimate of Box-Cox parameter based on Shapiro-Francia test statistic
out2$p.value # p.value of Shapiro-Francia test for transformed data 
out2$tf.data # transformed data set
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{functions}

