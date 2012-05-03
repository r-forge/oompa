\name{Mosaic}
\alias{Mosaic}
\alias{Mosaic-class}
\alias{plot,Mosaic,missing-method}
\alias{pltree,Mosaic-method}
\alias{summary,Mosaic-method}
\docType{class}
\title{Class "Mosaic"}
\description{
  Produce \dQuote{Eisen} plots of microarray or proteomics data.
}
\usage{
Mosaic(data,
       sampleMetric="pearson",
       sampleLinkage="average",
       geneMetric="euclid",
       geneLinkage="average",
       usecor=FALSE,
       center=FALSE,
       name="My mosaic")
\S4method{pltree}{Mosaic}(x, colors, labels, \dots)
\S4method{plot}{Mosaic,missing}(x, main=x@name, center=FALSE,
  scale=c("none", "row", "column"), limits=NULL,
  sampleColors=NULL, sampleClasses=NULL,
  geneColors=NULL, geneClasses=NULL, \dots)
}
\arguments{
  \item{data}{ Either a data frame or matrix with numeric values or an
    \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} as defined
    in the BioConductor tools for analyzing microarray data.}
  \item{sampleMetric}{Any valid distance metric that can be passed to the
    \code{\link{distanceMatrix}} function}
  \item{sampleLinkage}{Any valid linkage method that can be passed to the
    \code{\link{hclust}} function}
  \item{geneMetric}{Any valid distance metric that can be passed to the
    \code{\link{distanceMatrix}} function}
  \item{geneLinkage}{Any valid linkage method that can be passed to the
    \code{\link{hclust}} function}
  \item{center}{logical scalar. If \code{TRUE}, center the data rows.}
  \item{usecor}{logical scalar. If \code{TRUE}, scale the data rows to have
    standard deviation one.}
  \item{name}{character string specifying the name of this object}
  \item{x}{object of class \code{Mosaic}}
  \item{scale}{Same as in \code{\link{heatmap}}.}
  \item{colors}{An optional vector of character strings containing color
    names to be used when labeling the trees in the dendrogram. If
    provided, then the length should equal the number of columns in the
    original data matrix.}
  \item{labels}{An optional vector of character strings used to label
    the leaves in the dendrogram. If omitted, the column names are
    used.}
  \item{main}{character string specifying the main title for the plot}
  \item{limits}{An numeric vector. If provided, the data is truncated
    for display purposes, both above and below, at the minimum and
    maximum values of \code{limits}.}
  \item{sampleColors}{An optional character vector containing colors
    that will be used to label different sample types with a color bar
    across the top of the heat map.}
  \item{sampleClasses}{A logical vector or factor used to classify the
    samples into groups. Alternatively, an integer specifying the number
    \code{k} of groups into which to cut the sample dendrogram.}
  \item{geneColors}{An optional character vector containing colors
    that will be used to label different gene types with a color bar
    along the side of the heat map.}
  \item{geneClasses}{A logical vector or factor used to classify the
    genes into groups. Alternatively, an integer specifying the number
    \code{k} of groups into which to cut the gene dendrogram.}
  \item{\dots}{Additional parameters for \code{heatmap}.}
}
\details{
  One of the earliest papers in the microarray literature used
  independent clustering of the genes (rows) and samples (columns) to
  produce dendrograms that were plotted along with a red-green heat map
  of the centered expression values. Since that time, literally thousand
  of additional papers have published variations on these red-green
  images. R includes a function, \code{\link{heatmap}} that builds such
  figures.  However, that function is general purpose and has numerous
  optional parameters to tweak the display. The purpose of the
  \code{Mosaic} class is to provide a simplified object-oriented wrapper
  around \code{\link{heatmap}}, which as a side benefit allows us to
  keep track of the distance metrics and linkage rules that were used to
  produce the resulting figure.
}
\value{
  The \code{Mosaic} function constructs and returns a valid object of
  the \code{Mosaic} class.
}
\section{Objects from the Class}{
  Objects should be created with the \code{Mosaic} function.
}
\section{Slots}{
  \describe{
    \item{\code{data}:}{The \code{matrix} containing the numerical data }
    \item{\code{samples}:}{A dendrogram of class \code{hclust} produced
      by clustering the biological samples (columns of \code{data}). }
    \item{\code{genes}:}{A dendrogram of class \code{hclust} produced by
    clustering the genes (columns of \code{data}).}
    \item{\code{sampleMetric}:}{A \code{character} string; the distance
      metric used to cluster the samples. }
    \item{\code{sampleLinkage}:}{A \code{character} string; the linkage
      rule used to cluster the samples. }
    \item{\code{geneMetric}:}{A \code{character} string; the distance
      metric used to cluster the genes. }
    \item{\code{geneLinkage}:}{A \code{character} string; the linkage
      rule used to cluster the genes. }
    \item{\code{call}:}{An object of class \code{call} recording how the
    object was constructed. }
    \item{\code{name}:}{A \code{character} string; the name of this object. }
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = Mosaic, y = missing)}: Produce the
      \dQuote{Eisen} plot, using \code{\link{heatmap}}. }
    \item{pltree}{\code{signature(x = Mosaic)}: Plot the sample class
      dendrogram in the object. }
    \item{summary}{\code{signature(object = Mosaic)}: Write out a
      summary of the object. }
  }
}
\references{
  Eisen MB, Spellman PT, Brown PO, Botstein D.\cr
  "Cluster analysis and display of genome-wide expression patterns."\cr
  Proc Natl Acad Sci U S A. 1998 Dec 8;95(25):14863-8.
}
\author{
  Kevin R. Coombes \email{kcoombes@mdanderson.org},
  P. Roebuck \email{proebuck@mdanderson.org}
}
\seealso{
  \code{\link{cutree}},
  \code{\link{hclust}},
  \code{\link{heatmap}}
}
\examples{
showClass("Mosaic")

## simulate data from three different sample groups
d1 <- matrix(rnorm(100*10, rnorm(100, 0.5)), nrow=100, ncol=10, byrow=FALSE)
d2 <- matrix(rnorm(100*10, rnorm(100, 0.5)), nrow=100, ncol=10, byrow=FALSE)
d3 <- matrix(rnorm(100*10, rnorm(100, 0.5)), nrow=100, ncol=10, byrow=FALSE)
dd <- cbind(d1, d2, d3)
kind <- factor(rep(c('red', 'green', 'blue'), each=10))

## prepare the Mosaic object
m <- Mosaic(dd,
            sampleMetric='pearson',
            geneMetric='spearman',
            center=TRUE,
            usecor=TRUE)
summary(m)

## The default plot with red-green color map
plot(m, col=redgreen(64))

## change to a blue-yellow color map, and mark the four top splits in the
## sample direction with a color bar along the top
plot(m, col=blueyellow(128), sampleClasses=4,
     sampleColors=c('red', 'green', 'blue', 'black'))

## This time, mark the three classes that we know are there
plot(m, col=blueyellow(128), sampleClasses=kind,
     sampleColors=c('red', 'green', 'blue'))

plot(m, col=blueyellow(128),
     geneClasses=3, geneColors=c('red', 'green', 'black'))

## In addition, mark the top 5 splits in the gene dendrogram
plot(m,
     col=blueyellow(128),
     sampleClasses=kind,
     sampleColors=c('red', 'green', 'black'),
     geneClasses=5,
     geneColors=c('cyan', 'magenta', 'royalblue', 'darkgreen', 'orange'))

## plot the sample dendrogram by itself
cols <- as.character(kind)
pltree(m, labels=1:30, colors=cols)

## cleanup
rm(d1, d2, d3, dd, kind, cols, m)
}
\keyword{classes}
\keyword{cluster}
\keyword{hplot}
\keyword{multivariate}

