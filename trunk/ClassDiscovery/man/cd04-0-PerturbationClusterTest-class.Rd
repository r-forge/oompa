\name{PerturbationClusterTest}
\alias{PerturbationClusterTest}
\alias{PerturbationClusterTest-class}
\alias{summary,PerturbationClusterTest-method}
\docType{class}
\title{The PerturbationClusterTest Class }
\description{
  Performs a parametric bootstrap test (by adding independent Gaussian
  noise) to determine whether the clusters found by an unsupervised
  method appear to be robust in a given data set.
}
\usage{
PerturbationClusterTest(data, FUN, nTimes=100, noise=1, verbose=TRUE, \dots)
}
\arguments{
  \item{data}{A data matrix, numerical data frame, or
    \code{\link[Biobase]{ExpressionSet}} object.}
  \item{FUN}{A \code{function} that, given a data matrix,
    returns a vector of cluster assignments.  Examples of functions
    with this behavior are \code{\link{cutHclust}},
    \code{\link{cutKmeans}}, \code{\link{cutPam}}, and
    \code{\link{cutRepeatedKmeans}}.}
  \item{\dots}{Additional arguments passed to the classifying function, \code{FUN}.} 
  \item{noise}{An optional numeric argument; the standard deviation of
    the mean zero Gaussian noise added to each measurement during each
    bootstrap.  Defaults to 1.} 
  \item{nTimes}{The number of bootstrap samples to collect.}
  \item{verbose}{A logical flag}
}
\section{Objects from the Class}{
  Objects should be created using the  \code{PerturbationClusterTest}
  function, which performs the requested bootstrap on the
  clusters. Following the standard R paradigm, the resulting object can be
  summarized and plotted to determine the results of the test.
}
\section{Slots}{
  \describe{
    \item{\code{f}:}{A \code{function} that, given a data matrix,
      returns a vector of cluster assignments.  Examples of functions
      with this behavior are \code{\link{cutHclust}},
      \code{\link{cutKmeans}}, \code{\link{cutPam}}, and
      \code{\link{cutRepeatedKmeans}}. 
    }
    \item{\code{noise}:}{ The standard deviation of the Gaussian noise
      added during each bootstrap sample.}
    \item{\code{nTimes}:}{An integer, the number of bootstrap samples
      that were collected.}
    \item{\code{call}:}{An object of class \code{call}, which records
      how the object was produced. }
    \item{\code{result}:}{Object of class \code{matrix} containing, for
      each pair of columns in the original data, the number of times
      they belonged to the same cluster of a bootstrap sample. }
  }
}
\section{Extends}{
  Class \code{\link{ClusterTest}}, directly. See that class for
  descriptions of the inherited methods \code{image} and \code{hist}. 
}
\section{Methods}{
  \describe{
    \item{summary}{\code{signature(object = PerturbationClusterTest)}:
      Write out a summary of the object.}
  }
}
\references{
  Kerr MK, Churchill GJ.\cr
  \emph{Bootstrapping cluster analysis: Assessing the reliability of
  conclusions from microarray experiments.}\cr
  PNAS 2001; 98:8961-8965.
}
\author{
  Kevin R. Coombes \email{krc@silicovore.com}
}
\seealso{
  \code{\link{BootstrapClusterTest}},
  \code{\link{ClusterTest}}
}
\examples{
showClass("PerturbationClusterTest")

## simulate data from two different groups
d1 <- matrix(rnorm(100*30, rnorm(100, 0.5)), nrow=100, ncol=30, byrow=FALSE)
d2 <- matrix(rnorm(100*20, rnorm(100, 0.5)), nrow=100, ncol=20, byrow=FALSE)
dd <- cbind(d1, d2)
cols <- rep(c('red', 'green'), times=c(30,20))
colnames(dd) <- paste(cols, c(1:30, 1:20), sep='')
## peform your basic hierarchical clustering...
hc <- hclust(distanceMatrix(dd, 'pearson'), method='complete')

## bootstrap the clusters arising from hclust
bc <- PerturbationClusterTest(dd, cutHclust, nTimes=200, k=3, metric='pearson')
summary(bc)

## look at the distribution of agreement scores
hist(bc, breaks=101)

## let heatmap compute a new dendrogram from the agreement
image(bc, col=blueyellow(64), RowSideColors=cols, ColSideColors=cols)

## plot the agreement matrix with the original dendrogram
image(bc, dendrogram=hc, col=blueyellow(64), RowSideColors=cols, ColSideColors=cols)

## bootstrap the results of K-means
kmc <- PerturbationClusterTest(dd, cutKmeans, nTimes=200, k=3)
image(kmc, dendrogram=hc, col=blueyellow(64), RowSideColors=cols, ColSideColors=cols)

## contrast the behavior when all the data comes from the same group
xx <- matrix(rnorm(100*50, rnorm(100, 0.5)), nrow=100, ncol=50, byrow=FALSE)
hct <- hclust(distanceMatrix(xx, 'pearson'), method='complete')
bct <- PerturbationClusterTest(xx, cutHclust, nTimes=200, k=4, metric='pearson')
summary(bct)
image(bct, dendrogram=hct, col=blueyellow(64), RowSideColors=cols, ColSideColors=cols)

## cleanup
rm(d1, d2, dd, cols, hc, bc, kmc, xx, hct, bct)
}
\keyword{classes}
\keyword{multivariate}
\keyword{cluster}
\keyword{htest}

