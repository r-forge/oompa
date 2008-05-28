\name{matrixMean}
\alias{matrixMean}
\alias{matrixVar}
\title{Compute row-wise mean and variance using matrix operations}
\description{
  In large data sets, like those arising from microarray or proteomics
  experiments, it is often necessary to compute the mean and variance
  for each row among many thousands. The computations can be
  significantly faster if vectorized instead of using the inherent loop
  in an \code{apply} statement.
}
\usage{
matrixMean(x)
matrixVar(x, xmean)
}
\arguments{
  \item{x}{A matrix or data.frame containing numeric values.}
  \item{xmean}{A vector containing the means of the rows in \code{x}.}
}
\value{
  A vector containing a number of entries equal to the number of rows of
  \code{x}, where each entry is either the mean or the variance of the
  corresponding row.
}
\author{Kevin R. Coombes <kcoombes@mdanderson.org>}
\seealso{\code{\link{apply}}}
\examples{
nr <- 10000
nc <- 40
fake.data <- matrix(rnorm(nr*nc), nrow=nr)
fake.class <- rep(c('H', 'C'), each=20)

H.n <- sum(fake.class=='H')
C.n <- sum(fake.class=='C')

H.mean <- matrixMean(fake.data[, fake.class=='H'])
H.var <- matrixVar(fake.data[, fake.class=='H'], H.mean)

C.mean <- matrixMean(fake.data[, fake.class=='C'])
C.var <- matrixVar(fake.data[, fake.class=='C'], C.mean)

pooled.sd <- sqrt( (H.var*(H.n - 1) + C.var*(C.n - 1))/(H.n + C.n - 2) ) 

t.statistics <- (C.mean - H.mean)/pooled.sd/sqrt(1/H.n + 1/C.n)
}
\keyword{univar}
\keyword{array}
\keyword{algebra}

