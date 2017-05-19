\name{grades}
\alias{grades}
\docType{data}
\title{
Student Grades Data
}
\description{
Overall student grades for a class thaught by Dr. Ozlem Ilk
}
\usage{data(grades)}
\format{
  A data frame with 42 observations on the following variable.
  \describe{
    \item{\code{grades}}{a numeric vector for the student grades}
  }
}
%%  \details{
%%  ~~ If necessary, more details than the __description__ above %%  ~~
%%  }
%%  \source{
%%  ~~ reference to a publication or URL from which the data were %%  obtained ~~
%%  }
%%  \references{
%%  ~~ possibly secondary sources and usages ~~
%%  }
\examples{

library(AID)

data(grades)
hist(grades[,1])
out <- boxcoxnc(grades[,1])
confInt(out, level = 0.95)

}
\keyword{datasets}
