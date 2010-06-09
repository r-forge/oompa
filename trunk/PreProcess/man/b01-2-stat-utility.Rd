\name{stat.utility}
\alias{f.above.thresh}
\alias{f.cord}
\alias{f.oneway.rankings}

\title{OOMPA Statistical Utility Functions}
\description{
  Utility functions for statistical computations.
}
\usage{
f.above.thresh(a, t)
f.cord(x, y, inf.rm)
f.oneway.rankings(r, s)
}
\arguments{
  \item{a}{a vector}
  \item{t}{a real number}
  \item{x}{a vector}
  \item{y}{a vector}
  \item{inf.rm}{a logical value}
  \item{r}{vector}
  \item{s}{vector}
}
\value{
  \code{f.above.thresh} returns the fraction of elements in the vector
  \code{a} that are greater than the threshold \code{t}.

  \code{f.cord} returns the concordance coefficient between the two
  input vectors \code{x} and \code{y}. If \code{inf.rm} is true, then
  infinite values are removed before computing the concordance; missing
  values are always removed.

  \code{f.oneway.rankings} is implemented as \code{order(s)[r]} and I
  cannot recall why we defined it or where we used it.
}
\author{Kevin R. Coombes <kcoombes@mdanderson.org>}
\examples{
x <- rnorm(1000, 1, 2)
y <- rnorm(1000, 1, 2)
f.above.thresh(x, 0)
f.above.thresh(y, 0)
f.cord(x, y)
}
\keyword{univar}

