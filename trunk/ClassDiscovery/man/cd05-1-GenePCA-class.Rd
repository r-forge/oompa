\name{GenePCA}
\alias{GenePCA}
\alias{GenePCA-class}
\alias{plot,GenePCA,missing-method}
\docType{class}
\title{Class "GenePCA"}
\description{
  Perform principal components analysis on the genes (rows) from a
  microarray or proteomics experiment.
}
\usage{
GenePCA(geneData)
\S4method{plot}{GenePCA,missing}(x, splitter=0)
}
\arguments{
  \item{geneData}{A data matrix, with rows interpreted as genes and
    columns as samples}
  \item{x}{a \code{GenePCA} object}
  \item{splitter}{A logical vector classifying the genes.}
}
\details{
  This is a preliminary attempt at a class for principal components
  analysis of genes, parallel to the \code{\link{SamplePCA}} class for
  samples.  The interface will (one hopes) improve markedly in the next
  version of the library.
}
\value{
  The \code{GenePCA} function constructs and returns a valid object of
  the \code{GenePCA} class.
}
\section{Objects from the Class}{
  Objects should be created using the \code{GenePCA} generator function.
}
\section{Slots}{
  \describe{
    \item{\code{scores}:}{A \code{matrix} of size PxN, where P is the
      number of rows and N the number fo columns in the input,
      representing the projections of the input rows onto the first N
      principal components. }
    \item{\code{variances}:}{A \code{numeric} vector of length N; the
      amount of the total variance explained by each principal component.}
    \item{\code{components}:}{A \code{matrix} of size NxN containing
      each of the first P principal components as columns.}
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = GenePCA, y = missing)}: Plot the
      genes in the space of the first two principal components. }
  }
}
\author{
  Kevin R. Coombes \email{kcoombes@mdanderson.org}
}
\seealso{
  \code{\link{SamplePCA}},
  \code{\link{princomp}}
}
\examples{
# simulate samples from three different groups, with generic genes
d1 <- matrix(rnorm(100*10, rnorm(100, 0.5)), nrow=100, ncol=10, byrow=FALSE)
d2 <- matrix(rnorm(100*10, rnorm(100, 0.5)), nrow=100, ncol=10, byrow=FALSE)
d3 <- matrix(rnorm(100*10, rnorm(100, 0.5)), nrow=100, ncol=10, byrow=FALSE)
dd <- cbind(d1, d2, d3)

# perform PCA in gene space
gpc <- GenePCA(dd)

# plot the results
plot(gpc)

# cleanup
rm(d1, d2, d3, dd, gpc)
}
\keyword{algebra}
\keyword{array}
\keyword{classes}
\keyword{multivariate}

