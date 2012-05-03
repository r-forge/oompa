\name{distanceMatrix}
\alias{distanceMatrix}
\title{Distance Matrix Computation}
\description{
  This function computes and returns the distance matrix determined by
  using the specified distance metric to compute the distances between
  the columns of a data matrix.
}
\usage{
distanceMatrix(dataset, metric, \dots)
}
\arguments{
  \item{dataset}{A numeric matrix or an \code{\link[Biobase]{ExpressionSet}}}
  \item{metric}{ A character string defining the distance metric. This
    can be \code{pearson}, \code{sqrt pearson}, \code{spearman},
    \code{absolute pearson}, \code{uncentered correlation}, \code{weird}
    or any of the metrics accepted by the     \code{\link{dist}}
    function. At present, the latter function accepts
    \code{euclidean}, \code{maximum}, \code{manhattan}, \code{canberra},
    \code{binary}, or \code{minkowski}. Any initial substring that
    uniquely defines one of the metrics will work.}
  \item{\dots}{Additional parameters to be passed on to
    \code{\link{dist}}. }
}
\details{
  This function differs from \code{\link{dist}} in two ways, both of
  which are motivated by common practice in the analysis of microarray
  or proteomics data. First, it computes distances between column vectors
  instead of between row vectors. In a typical microarray experiment,
  the data is organized so the rows represent genes and the columns
  represent different biological samples. In many applications,
  relations between the biological samples are more interesting than
  relationships between genes. Second, \code{distanceMatrix} adds
  additional distance metrics based on correlation.
  \itemize{
    \item{\code{pearson}}{The most common metric used in the microarray literature is
      the \code{pearson} distance, which can be computed in terms of the
      Pearson correlation coefficient as \code{(1-cor(dataset))/2}.}
    \item{\code{uncentered correlation}}{This metric was introduced in
      the Cluster and TreeView software from the Eisen lab at
      Stanford. It is computed using the formulas for Pearson
      correlation, but assuming that both vectors have mean zero.}
    \item{\code{spearman}}{The \code{spearman} metric used the same formula, but
      substitutes the Spearman rank correlation for the Pearson
      correlation.}
    \item{\code{absolute pearson}}{The \code{absolute pearson} metric used the absolute
      correlation coefficient; i.e., \code{(1-abs(cor(dataset)))}.}
    \item{\code{sqrt pearson}}{The \code{sqrt pearson} metric used the square root of the
      pearson distance metric; i.e., \code{sqrt(1-cor(dataset))}.}
    \item{\code{weird}}{The \code{weird} metric uses the Euclidean distance between
      the vectors of correlation coefficients; i.e., dist(cor(dataset)).}
  }
}
\value{
  A distance matrix in the form of an object of class \code{dist}, of
  the sort returned by the \code{dist} function or the \code{as.dist}
  function. 
}
\author{
  Kevin R. Coombes \email{kcoombes@mdanderson.org}
}
\seealso{
  \code{\link{dist}},
  \code{\link{as.dist}}
}
\section{BUGS}{
  It would be good to have a better name for the \code{weird} metric.
}
\examples{
dd <- matrix(rnorm(100*5, rnorm(100)), nrow=100, ncol=5)
distanceMatrix(dd, 'pearson')
distanceMatrix(dd, 'euclid')
distanceMatrix(dd, 'sqrt')
distanceMatrix(dd, 'weird')
rm(dd) # cleanup
}
\keyword{array}

