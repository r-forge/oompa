\name{justClusters}
\alias{cutHclust}
\alias{cutPam}
\alias{cutKmeans}
\alias{repeatedKmeans}
\alias{cutRepeatedKmeans}

\title{Get the List of Classes From A Clustering Algorithm }
\description{
  Unsupervised clustering algorithms, such as partitioning around medoids
  (\code{\link[cluster]{pam}}), K-means (\code{\link{kmeans}}), or
  hierarchical clustering (\code{\link{hclust}}) after cutting the tree,
  produce a list of class assignments along with other structure. To
  simplify the interface for the \code{\link{BootstrapClusterTest}}  and
  \code{\link{PerturbationClusterTest}}, we have written these routines
  that simply extract these cluster assignments.
}
\usage{
cutHclust(data, k, method = "average", metric = "pearson")
cutPam(data, k)
cutKmeans(data, k)
cutRepeatedKmeans(data, k, nTimes)

repeatedKmeans(data, k, nTimes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A numerical data matrix}
  \item{k}{The number of classes desired from the algorithm}
  \item{method}{Any valid linkage method that can be passed to the
    \code{\link{hclust}} function}
  \item{metric}{Any valid distance metric that can be passed to the
    \code{\link{distanceMatrix}} function}
  \item{nTimes}{An integer; the number of times to repeat the K-means
    algorithm with a different random starting point}
}
\details{
  Each of the clustering routines used here has a different
  structure for storing cluster assignments. The \code{\link{kmeans}}
  function stores the assignments in a \sQuote{cluster} attribute.  The
  \code{\link[cluster]{pam}} function uses a \sQuote{clustering} attribute.  For
  \code{\link{hclust}}, the assigments are produced by a call to the
  \code{\link{cutree}} function.

  It has been observed that the K-means algorithm can converge to
  different solutions depending on the starting values of the group
  centers. We also include a routine (\code{repeatedKmeans}) that runs
  the K-means algorithm repeatedly, using different randomly generated
  staring points each time, saving the best results.
}
\value{
  Each of the \code{cut...} functions returns a vector of integer values
  representing the cluster assignments found by the algorithm.

  The \code{repeatedKmeans} function returns a list \code{x} with three
  components.  The component \code{x$kmeans} is the result of the call
  to the \code{kmeans} function that produced the best fit to the
  data. The component \code{x$centers} is a matrix containing the list
  of group centers that were used in the best call to \code{kmeans}.
  The component \code{x$withinss} contains the sum of the within-group
  sums of squares, which is used as the measure of fitness.
}
\author{Kevin R. Coombes <kcoombes@mdanderson.org>}
\seealso{
  \code{\link[cluster]{pam}}, \code{\link{hclust}},
  \code{\link{cutree}}, \code{\link{kmeans}}. 
}
\examples{
# simulate data from three different groups
d1 <- matrix(rnorm(100*10, rnorm(100, 0.5)), nrow=100, ncol=10, byrow=FALSE)
d2 <- matrix(rnorm(100*10, rnorm(100, 0.5)), nrow=100, ncol=10, byrow=FALSE)
d3 <- matrix(rnorm(100*10, rnorm(100, 0.5)), nrow=100, ncol=10, byrow=FALSE)
dd <- cbind(d1, d2, d3)

cutKmeans(dd, k=3)
cutKmeans(dd, k=4)

cutHclust(dd, k=3)
cutHclust(dd, k=4)

cutPam(dd, k=3)
cutPam(dd, k=4)

cutRepeatedKmeans(dd, k=3, nTimes=10)
cutRepeatedKmeans(dd, k=4, nTimes=10)

# cleanup
rm(d1, d2, d3, dd)
}
\keyword{cluster}
\keyword{multivariate}
