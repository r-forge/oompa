\name{matrixT}
\alias{matrixT}
\alias{matrixUnequalT}
\alias{matrixPairedT}
\alias{matrixVar}
\alias{matrixMean}

\title{ Multiple T Tests by Matrix Multiplication }
\description{
  Utility functions for computing vectors of row-by-row means,
  variances, and t-statistics.
}
\usage{
matrixMean(x,na.rm=FALSE)
matrixVar(x, xmean, na.rm=FALSE)
matrixT(m, v, na.rm=FALSE)
matrixUnequalT(m, v)
matrixPairedT(m, v, pf)
}
\arguments{
  \item{x}{ a matrix }
  \item{xmean}{ a numeric vector or single-column matrix }
  \item{m}{ a matrix }
  \item{na.rm}{ a logical value indicating whether means, variances, and
    t-statistics should be computed after omitting NA values
    from individual rows of the data matrix.}
  \item{v}{ a logical vector of length equal to the number of columns of \code{m} }
  \item{pf}{ a numerical vector of length equal to the number of columns
    of \code{m}, indicating which samples should be paired }
}
\value{
  \code{matrixMean} returns a single-column matrix containing the
  row-by-row means of \code{x}.

  \code{matrixVar} returns a single-column matrix containing the
  row-by-row means of \code{x}, assuming that \code{xmean} contains the
  corresponding mean values.

  \code{matrixT} returns a single-column matrix of t-statistics from a
  two-sample t-test comparing the columns for which \code{v} is true to
  those for which \code{v} is false.

  \code{matrixPairedT} returns a single-column matrix of t-statistics from a
  paired t-test.

  \code{matrixUnequalT} returns a list with two components: \code{tt} is
  a single-column matrix of t-statistics from a two-sample unequal
  variance t-test comparing the columns for which \code{v} is true to
  those for which \code{v} is false, and \code{df} is a single-column
  matrix of the degrees of freedom associated with each row.. 
}

\author{Kevin R. Coombes <krc@silicovore.com>}
\examples{
ng <- 1000
ns <- 50
dat <- matrix(rnorm(ng*ns), ncol=ns)
clas <- factor(rep(c('Good', 'Bad'), each=25))
myMean <- matrixMean(dat)
myVar  <- matrixVar(dat, myMean)
plot(myMean, myVar)


myT <- matrixT(dat, clas)
hist(myT)

rm(ng, ns, dat, myMean, myVar, myT)
}
\keyword{array}
\keyword{htest}
\keyword{multivariate}
