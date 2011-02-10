\name{BootstrapClusterTest}
\docType{class}
\alias{BootstrapClusterTest}
\alias{BootstrapClusterTest-class}
\alias{summary,BootstrapClusterTest-method}

\title{The BootstrapClusterTest Class }
\description{
  Performs a nonparametric bootstrap (sampling with replacement) test to
  determine whether the clusters found by an unsupervised method appear
  to be robust in a given data set.
}
\usage{
BootstrapClusterTest(data, FUN, subsetSize, nTimes = 100, verbose = TRUE, ...)
}
\arguments{
  \item{data}{A data matrix, numerical data frame, or
    \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} object.}
  \item{FUN}{A \code{function} that, given a data matrix,
    returns a vector of cluster assignments.  Examples of functions
    with this behavior are \code{\link{cutHclust}},
    \code{\link{cutKmeans}}, \code{\link{cutPam}}, and
    \code{\link{cutRepeatedKmeans}}.}
  \item{\dots}{Additional arguments passed to the classifying function, \code{FUN}.} 
  \item{subsetSize}{An optional integer argument. If present,
      each iteration of the bootstrap selects \code{subsetSize} rows
      from the original data matrix. If missing, each bootstrap contains
    the same number of rows as the original data matrix.} 
  \item{nTimes}{The number of bootstrap samples to collect.}
  \item{verbose}{A logical flag}
}
\section{Objects from the Class}{
Objects should be created using the  \code{BootstrapClusterTest}
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
    \item{\code{subsetSize}:}{ The number of rows to be included in each
      bootstrap sample.}
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
    \item{summary}{\code{signature(object = BootstrapClusterTest)}:
      Write out a summary of the object. }
  }
}
\references{
  Kerr MK, Churchill GJ. Boostrapping cluster analysis: Assessing the
  reliability of conclusions from microarray experiments. PNAS 2001;
  98:8961-8965.
}
\author{Kevin R. Coombes <kcoombes@mdanderson.org>}
\seealso{
  \code{\link{ClusterTest}},
  \code{\link{PerturbationClusterTest}}
}
\examples{
# simulate data from two different groups
d1 <- matrix(rnorm(100*30, rnorm(100, 0.5)), nrow=100, ncol=30, byrow=FALSE)
d2 <- matrix(rnorm(100*20, rnorm(100, 0.5)), nrow=100, ncol=20, byrow=FALSE)
dd <- cbind(d1, d2)
cols <- rep(c('red', 'green'), times=c(30,20))
# peform your basic hierarchical clustering...
hc <- hclust(distanceMatrix(dd, 'pearson'), method='complete')

# bootstrap the clusters arising from hclust
bc <- BootstrapClusterTest(dd, cutHclust, nTimes=200, k=3, metric='pearson')
summary(bc)

# look at the distribution of agreement scores
hist(bc, breaks=101)

# let heatmap compute a new dendrogram from the agreement
image(bc, col=blueyellow(64), RowSideColors=cols, ColSideColors=cols)

# plot the agreement matrix with the original dendrogram
image(bc, dendrogram=hc, col=blueyellow(64), RowSideColors=cols, ColSideColors=cols)

# bootstrap the results of PAM
pamc <- BootstrapClusterTest(dd, cutPam, nTimes=200, k=3)
image(pamc, dendrogram=hc, col=blueyellow(64), RowSideColors=cols, ColSideColors=cols)

# contrast the behavior when all the data comes from the same group
xx <- matrix(rnorm(100*50, rnorm(100, 0.5)), nrow=100, ncol=50, byrow=FALSE)
hct <- hclust(distanceMatrix(xx, 'pearson'), method='complete')
bct <- BootstrapClusterTest(xx, cutHclust, nTimes=200, k=4, metric='pearson')
summary(bct)
image(bct, dendrogram=hct, col=blueyellow(64), RowSideColors=cols, ColSideColors=cols)

# cleanup
rm(d1, d2, dd, cols, hc, bc, pamc, xx, hct, bct)
}
\keyword{classes}
\keyword{multivariate}
\keyword{cluster}
\keyword{htest}
