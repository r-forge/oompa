\name{SamplePCA}
\alias{SamplePCA}
\alias{SamplePCA-class}
\alias{plot,SamplePCA,missing-method}
\alias{predict,SamplePCA-method}
\alias{screeplot,SamplePCA-method}
\alias{summary,SamplePCA-method}
\alias{identify,SamplePCA-method}
\alias{text,SamplePCA-method}
\docType{class}
\title{Class "SamplePCA"}
\description{
  Perform principal components analysis on the samples (columns) from a
  microarray or proteomics experiment.
}
\usage{
SamplePCA(data, splitter=0, usecor=FALSE, center=TRUE)
\S4method{plot}{SamplePCA,missing}(x, splitter=x@splitter, col, main='', which=1:2, \dots)
}
\arguments{
  \item{data}{ Either a data frame or matrix with numeric values or an
    \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} as defined
    in the BioConductor tools for analyzing microarray data.}
  \item{splitter}{ If \code{data} is a data frame or matrix, then splitter
    must be either a logical vector or a factor. If \code{data} is an
    \code{ExpressionSet}, then \code{splitter} can be a character string that
    names one of the factor columns in the associated
    \code{\link[Biobase]{phenoData}} subobject.}
  \item{center}{A logical value; should the rows of the data matrix be
    centered first?}
  \item{usecor}{A logical value; should the rows of the data matrix be
    scaled to have standard deviation 1?}
  \item{x}{A \code{SamplePCA} object}
  \item{col}{A list of colors to represent each level of the
    \code{splitter} in the plot. If this parameter is missing, the
    function will select colors automatically.}
  \item{main}{A character string; the plot title}
  \item{which}{A numeric vector of length two specifying which two principal
    components should be included in the plot.}
  \item{\dots}{Additional graphical parameters for \code{plot}}.
}
\section{Objects from the Class}{
  Objects should be created using the \code{SamplePCA} function. In the
  simplest case, you simply pass in a data matrix and a logical vector,
  \code{splitter}, assigning classes to the columns, and the constructor
  performs principal components analysis on the column.  The
  \code{splitter} is ignored by the constructor and is simply saved to
  be used by the plotting routines.  If you omit the \code{splitter},
  then no grouping structure is used in the plots.

  If you pass \code{splitter} as a factor instead of a logical vector,
  then the plotting routine will distinguish all levels of the factor.
  The code is likely to fail, however, if one of the levels of the
  factor has zero representatives among the data columns.

  As with the class comparison functions (see, for example,
  \code{\link[ClassComparison]{MultiTtest}}) that are part of OOMPA,
  we can also perform PCA on
  \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} objects 
  from the BioConductor libraries. In this case, we pass in an
  \code{ExpressionSet} object along with a character string containing the
  name of a factor to use for splitting the data.
}
\section{Slots}{
  \describe{
    \item{\code{scores}:}{A \code{matrix} of size NxN, where N is the
      number of columns in the input, representing the projections of
      the input columns onto the first N principal components. }
    \item{\code{variances}:}{A \code{numeric} vector of length N; the
      amount of the total variance explained by each principal component.}
    \item{\code{components}:}{A \code{matrix} of size PxN (the same size
      as the input matrix) containing each of the first P principal
      components as columns.}
    \item{\code{splitter}:}{A logical vector or factor of length N
      classifying the columns into known groups.}
    \item{\code{usecor}:}{A \code{logical} value; was the data standardized?}
    \item{\code{shift}:}{A \code{numeric} vector of length P; the mean
      vector of the input data, which is used for centering by the
      \code{predict} method.}
    \item{\code{scale}:}{A \code{numeric} vector of length P; the
      standard deviation of the input data, which is used for scaling by
    the \code{predict} method.}
    \item{\code{call}:}{An object of class \code{call} that records
      how the object was created. }
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = SamplePCA, y = missing)}: Plot the
      samples in a two-dimensional principal component space. }
    \item{predict}{\code{signature(object = SamplePCA)}: Project new
      data into the principal component space. }
    \item{screeplot}{\code{signature(x = SamplePCA)}: Produce a bar
      chart of the variances explained by each principal component. }
    \item{summary}{\code{signature(object = SamplePCA)}: Write out a
      summary of the object.}
    \item{identify}{\code{signature(object = SamplePCA)}: interactively
      identify points in the plot of a \code{SamplePCA} object.}
    \item{text}{\code{signature(object = SamplePCA)}: Add sample
      identifiers to the scatter plot of a \code{SamplePCA} object,
      using the base \code{text} function.}
  }
}
\details{
  The main reason for developing the \code{SamplePCA} class is that the
  \code{\link{princomp}} function is very inefficient when the number of
  variables (in the microarray setting, genes) far exceeds the number of
  observations (in the microarray setting, biological samples).  The
  \code{princomp} function begins by computing the full covariance
  matrix, which gets rather large in a study involving tens of thousands
  of genes.  The \code{SamplePCA} class, by contrast, uses singular
  value decomposition (\code{\link{svd}}) on the original data matrix to
  compute the principal components.

  The base functions \code{screeplot}, which produces  a barplot of the
  percentage of variance explained by each component, and \code{plot},
  which produces a scatter plot comparing two selected components
  (defaulting to the first two), have been generalized as methods for
  the \code{SamplePCA} class.  You can add sample labels to the scatter
  plot using either the \code{text} or \code{identify} methods.  One
  should, however, note that the current implementaiton of these methods
  only works when plotting the first two components.
}
\value{
  The \code{SamplePCA} function constructs and returns an object of the
  \code{SamplePCA} class.  We assume that the input data matrix has N
  columns (of biological samples) and P rows (of genes). 

  The \code{predict} method returns a matrix whose size is the number of
  columns in the input by the number of principal components.
}
\author{
  Kevin R. Coombes \email{krc@silicovore.com}
}
\seealso{
  \code{\link{princomp}},
  \code{\link{GenePCA}}
}
\examples{
showClass("SamplePCA")

## simulate data from three different groups
d1 <- matrix(rnorm(100*10, rnorm(100, 0.5)), nrow=100, ncol=10, byrow=FALSE)
d2 <- matrix(rnorm(100*10, rnorm(100, 0.5)), nrow=100, ncol=10, byrow=FALSE)
d3 <- matrix(rnorm(100*10, rnorm(100, 0.5)), nrow=100, ncol=10, byrow=FALSE)
dd <- cbind(d1, d2, d3)
kind <- factor(rep(c('red', 'green', 'blue'), each=10))
colnames(dd) <- paste(kind, rep(1:10, 3), sep='')

## perform PCA
spc <- SamplePCA(dd, splitter=kind)

## plot the results
plot(spc, col=levels(kind))

## mark the group centers
x1 <- predict(spc, matrix(apply(d1, 1, mean), ncol=1))
points(x1[1], x1[2], col='red', cex=2)
x2 <- predict(spc, matrix(apply(d2, 1, mean), ncol=1))
points(x2[1], x2[2], col='green', cex=2)
x3 <- predict(spc, matrix(apply(d3, 1, mean), ncol=1))
points(x3[1], x3[2], col='blue', cex=2)

## check out the variances
screeplot(spc)

## cleanup
rm(d1, d2, d3, dd,kind, spc, x1, x2, x3)
}
\keyword{classes}
\keyword{array}
\keyword{algebra}
\keyword{multivariate}
