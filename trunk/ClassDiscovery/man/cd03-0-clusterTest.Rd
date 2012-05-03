\name{ClusterTest-class}
\docType{class}
\alias{ClusterTest}
\alias{ClusterTest-class}
\alias{hist,ClusterTest-method}
\alias{image,ClusterTest-method}
\alias{summary,ClusterTest-method}

\title{The ClusterTest Class}
\description{ This is a base class for tests that attempt to determine
  whether the groups found by an unsupervised clustering method are
  statistically significant.
}
\usage{
\S4method{image}{ClusterTest}(x, dendrogram, ...)
}
\arguments{
  \item{x}{An object of the \code{ClusterTest} class.}
  \item{dendrogram}{An object with S3 class \code{hclust}, as returned
    by the \code{\link{hclust}} function.}
  \item{\dots}{Additional graphical parameters to be passed to the
    standard \code{heatmap} function that is used to implement the
    \code{image} method.}
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ClusterTest",  ...)}. 
Most users, however, will only create objects from one of the derived
classes such as \code{\link{BootstrapClusterTest}} or
\code{\link{PerturbationClusterTest}}. 
}
\section{Slots}{
  \describe{
    \item{\code{call}:}{An object of class \code{call}, which shows
      how the object was constructed.}
    \item{\code{result}:}{A symmetric \code{matrix} containing the
      results of the cluster reproducibility test.  The size of the
      matrix corresponds to the number of samples (columns) in the data
      set on which the test was performed.  The \code{result} matrix
      should contain "agreement" values between 0 and 1, representing for
      each pair of samples the fraction of times that they were
      collected into the same cluster.} 
  }
}
\section{Methods}{
  \describe{
    \item{hist}{\code{signature(x = "ClusterTest")}: Produces a
      histogram of the agreement fractions.  When a true group structure
      exists, one expects a multimodal distribution,with one group of
      agreements near 0 (for pairs belonging to different clusters) and
      one group of agreements near 1 (for pairs belonging to the same
      cluster).} 
    \item{image}{\code{signature(x = "ClusterTest")}: Uses the
      \code{heatmap} function to display the agreement matrix. The
      optional \code{dendrogram} argument should be used to display the
      extent to which the agreement matrix matches the results of
      hierarchical clustering using the full data set.  This method
      invisibly returns the result of a call to \code{heatmap}; thus, you
      can use \code{keep.dendro=TRUE} to recover the cluster assignments
      from the dendrograms.}
    \item{summary}{\code{signature(object = "ClusterTest")}: Write out a
    summary of the object.}
  }
}
\references{
  Kerr MK, Churchill GJ. Bootstrapping cluster analysis: Assessing the
  reliability of conclusions from microarray experiments. PNAS 2001;
  98:8961-8965.
}
\author{Kevin R. Coombes <kcoombes@mdanderson.org>}
\seealso{
  \code{\link{BootstrapClusterTest}},
  \code{\link{PerturbationClusterTest}}, \code{\link{heatmap}}
}
\examples{
# simulate data from two different classes
d1 <- matrix(rnorm(100*30, rnorm(100, 0.5)), nrow=100, ncol=30, byrow=FALSE)
d2 <- matrix(rnorm(100*20, rnorm(100, 0.5)), nrow=100, ncol=20, byrow=FALSE)
dd <- cbind(d1, d2)

# cluster the data
hc <- hclust(distanceMatrix(dd, 'pearson'), method='average')

# make a fake reproducibility matrix
fraud <- function(x) {
  new('ClusterTest', result=abs(cor(x)), call=match.call())
}

fake <- fraud(dd)
summary(fake)

hist(fake)

image(fake) # let heatmap compute a new dendrogram from the agreements

image(fake, dendrogram=hc) # use the actual dendrogram from the data

image(fake, dendrogram=hc, col=blueyellow(64)) # change the colors

#cleanup
rm(fake, fraud, hc, dd, d1, d2)
}
\keyword{classes}
\keyword{multivariate}
\keyword{cluster}
\keyword{htest}
