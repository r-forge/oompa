\name{cluster3}
\alias{cluster3}
\title{ Cluster a Dataset Three Ways }
\description{
  Produces and plots dendrograms using three similarity measures:
  Euclidean distance, Pearson correlation, and Manhattan distance on
  dichotomized data.
}
\usage{
cluster3(data, eps=logb(1, 2), name="", labels=dimnames(data)[[2]])
}
\arguments{
  \item{data}{ A matrix, numeric data.frame, or
    \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} object.}
  \item{eps}{ A numerical value; the threshold at which to dichotomize
    the data}
  \item{name}{A character string to label the plots }
  \item{labels}{A vector of character strings used to label the items in
    the dendrograms.}
}
\value{
  Invisibly returns the \code{data} object on which it was invoked.
}
\author{
  Kevin R. Coombes email{krc@silicovore.com}
}
\seealso{
  \code{\link{hclust}}
}
\examples{
## simulate data from two different classes
d1 <- matrix(rnorm(100*30, rnorm(100, 0.5)), nrow=100, ncol=30, byrow=FALSE)
d2 <- matrix(rnorm(100*20, rnorm(100, 0.5)), nrow=100, ncol=20, byrow=FALSE)
dd <- cbind(d1, d2)
## cluster it 3 ways
par(mfrow=c(2,2))
cluster3(dd)
par(mfrow=c(1,1))
}
\keyword{hplot}
\keyword{cluster}

