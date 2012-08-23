\name{boxcoxnc}
\alias{boxcoxnc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Function to Estimate Box-Cox Power Transformation Parameter via Normality Tests and Artifical Covariate Method
}
\description{
\code{boxcoxnc} utilizes seven different normality tests and artifical covariate method to estimate Box-Cox power transformation parameter and provides graphical analysis
}
\usage{boxcoxnc(data, method = "all", lam = seq(-2, 2, 0.01), plotit = TRUE, rep = 30)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{is a vector, matrix for univariate dataset}
  \item{method}{expects a character string to select the desired method to be used to estimate Box-Cox transformation parameter. To use Shapiro-Wilk test method should be set to "sw". For method = "ad", boxcoxnc function uses Anderson-Darling test to estimste Box-Cox transformation parameter. Similarly, method should be set to "cvm", "pt", "sf", "lt", "jb", "ac" to use Cramer-von Mises, Pearson Chi-square, Shapiro-Francia, Lilliefors, Jarque-Bera tests and artificial covariate method, respectively. To use all the methods at the same time, one should use method = "ad".}
  \item{lam}{is a vector which includes the sequence of candidate lambda values. Default is set to (-2,2) with increment 0.01} 
  \item{plotit}{plots normality test statistic vs lambda for methods utilizing normality tests. It also plots log-likelihood vs lambda for artificial covariate method. Defaults plotit = TRUE} 
  \item{rep}{is an integer which denotes the replication number for artificial covariate method. Default is set to 30}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Returns a matrix of output with the results of seven different normality tests and artifical covariate method. The first row of the matrix corresponds to the related estimates of lambda. The subsequent rows correspond to the p-values of different normality tests for each estimates of lambda.
}
\references{
Box, G. E. P., Cox, D. R. (1964). An Analysis of Transformations (with discussion). \emph{Journal of the Royal Statistical Society, Series B (Methodological)}, \bold{26}, 211--252.

Dag, O., Asar, O., Ilk, O. (2012). A Methodology to Implement Box-Cox Transformation When No Covariate is Available. Submitted.

Gross, J. (2012). nortest: Tests for Normality. R package version 1.0.

Thode, H. C. (2002). \emph{Testing for Normality}. New York: Marcel Dekker.

Trapletti, A., Hornik, K. (2012). tseries: Time Series Analysis and Computational Finance. R package version. 0.10-28.

Venables, W. N., Ripley, B. D. (2002). \emph{Modern Applied Statistics with S 4th ed.}. New York: Springer.

}
\author{
Osman Dag, Ozgur Asar, Ozlem Ilk
}
\note{
This is the version 1.0 of this user documentation file.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%\code{\link{boxcox}}
%\code{\link{powerTransform}}
%}
\examples{
data(textile)
boxcoxnc(textile[,1])
boxcoxnc(textile[,1])$result[1,1]
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ power transformation }

