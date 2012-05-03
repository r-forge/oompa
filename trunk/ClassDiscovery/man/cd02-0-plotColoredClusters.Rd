\name{plotColoredClusters}
\alias{plotColoredClusters}
\alias{pcc}
\title{Plot Dendrograms with Color-Coded Labels}
\description{
  Provides an interface to \code{\link{plclust}} that makes it easier to
  plot dendrograms with labels that are color-coded, usually to indicate
  the different levels of a factor.
}
\usage{
plotColoredClusters(hd, labs, cols, cex = 0.8, main = "", line = 0, ...)
pcc(hd, colors=NULL, ...)
}
\arguments{
  \item{hd}{ An object with S3 class \code{hclust}, as produced by the
    \code{hclust} function.}
  \item{labs}{ A vector of character strings used to label the leaves
    in the dendrogram}.
  \item{cols}{ A vector of color names suitable for passing to the
    \code{col} argument of graphics routines.}
  \item{cex}{ A numeric value; the character expansion parameter of
    \code{par}.}
  \item{main}{ A character string; the plot title}
  \item{line}{ An integer determining how far away to plot the labels;
    see \code{\link{mtext}} for details.}
  \item{colors}{ A list; see details.}
  \item{\dots}{Any additional graphical parameters that can be supplied
    when plotting an \code{hclust} object.}
}
\details{
  The \code{plotColoredClusters} function is used to implement the
  \code{pltree} methods of the 
  \code{Mosaic} class and the \code{PCanova} class.  It simply bundles
  a two step process (first plotting the dendrogram with no labels,
  followed by writing the labels in the right places with the desired
  colors) into a single unit.

  The \code{pcc} function also produces dendrograms with colored
  annotations.  However, instead of coloring the labels based on a
  single factor, it produces color bars for any number of factors.  The
  \code{colors} argument should be a list with named components, where
  each component should correspond to a factor and a color scheme.
  Specifically, the components must themselves be lists with two
  components named \code{fac} (and containing the factor) and \code{col}
  (containing a named vector specifying colors for each level of the
  factor).
}
\value{
  The function has no useful return value; it merely produces a plot.
}
\seealso{ \code{\link{plclust}}, \code{\link{hclust}},
  \code{\link{Mosaic}}, \code{\link{PCanova}}, \code{\link{par}}
}
\examples{
# simulate data from three different groups
d1 <- matrix(rnorm(100*10, rnorm(100, 0.5)), nrow=100, ncol=10, byrow=FALSE)
d2 <- matrix(rnorm(100*10, rnorm(100, 0.5)), nrow=100, ncol=10, byrow=FALSE)
d3 <- matrix(rnorm(100*10, rnorm(100, 0.5)), nrow=100, ncol=10, byrow=FALSE)
dd <- cbind(d1, d2, d3)

# perform hierarchical clustering using correlation
hc <- hclust(distanceMatrix(dd, 'pearson'), method='average')
cols <- rep(c('red', 'green', 'blue'), each=10)
labs <- paste('X', 1:30, sep='')

# plot the dendrogram with color-coded groups
plotColoredClusters(hc, labs=labs, cols=cols)

# simulate another dataset
fakedata <- matrix(rnorm(200*30), ncol=30)
colnames(fakedata) <- paste("P", 1:30, sep='')
# define two basic factors, with colors
faccol <- list(fac=factor(rep(c("A", "B"), each=15)),
               col=c(A='red', B='green'))
fac2col <- list(fac=factor(rep(c("X", "Y", "Z"), times=10)),
               col=c(X='cyan', Y='orange', Z='purple'))
# add another factor that reverses the colors
BA <- faccol
BA$col <- c(A='blue', B='yellow')
# assemble the list of factors
colors <- list(AB=faccol, XYZ=fac2col, "tricky long name"=fac2col,
               another=BA)
# cluster the samples
hc <- hclust(distanceMatrix(fakedata, "pearson"), "ward")
# plot the results
pcc(hc, colors)


#cleanup
rm(d1, d2, d3, dd, hc, cols, labs)
rm(fakedata, faccol, fac2col, BA, colors)
}
\keyword{hplot}
\keyword{cluster}
