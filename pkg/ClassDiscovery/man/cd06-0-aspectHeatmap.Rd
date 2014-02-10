\name{aspectHeatmap}
\alias{aspectHeatmap}
\title{Heatmap with control over the aspect ratio}
\description{
  This function replaces the heatmap function in the stats package with
  a function with additional features. First, the user can specify an
  aspect ratio instead of being forced to accept a square image even
  when the matrix is not square. Second, the user can overlay horizontal
  or vertical lines to mark out important regions.
}
\usage{
aspectHeatmap(x,
              Rowv=NULL,
              Colv=if (symm) "Rowv" else NULL,
              distfun=dist,
              hclustfun=hclust,
              reorderfun=function(d, w) reorder(d, w),
              add.expr,
              symm=FALSE,
              revC=identical(Colv, "Rowv"),
              scale=c("row", "column", "none"),
              na.rm=TRUE,
              margins=c(5, 5),
              ColSideColors,
              RowSideColors,
              cexRow=0.2 + 1/log10(nr),
              cexCol=0.2 + 1/log10(nc),
              labRow=NULL,
              labCol=NULL,
              main=NULL,
              xlab=NULL,
              ylab=NULL,
              keep.dendro=FALSE,
              verbose=getOption("verbose"),
              hExp=1,
              wExp=1,
              vlines=NULL,
              hlines=NULL,
              lasCol=2,
              \dots)
}
\arguments{
  \item{x}{ See documentation for \link{heatmap}.}
  \item{Rowv}{ See documentation for \link{heatmap}.}
  \item{Colv}{ See documentation for \link{heatmap}.}
  \item{distfun}{ See documentation for \link{heatmap}.}
  \item{hclustfun}{ See documentation for \link{heatmap}.}
  \item{reorderfun}{ See documentation for \link{heatmap}.}
  \item{add.expr}{ See documentation for \link{heatmap}.}
  \item{symm}{ See documentation for \link{heatmap}.}
  \item{revC}{ See documentation for \link{heatmap}.}
  \item{scale}{ See documentation for \link{heatmap}.}
  \item{na.rm}{ See documentation for \link{heatmap}.}
  \item{margins}{ See documentation for \link{heatmap}.}
  \item{ColSideColors}{ See documentation for \link{heatmap}.}
  \item{RowSideColors}{ See documentation for \link{heatmap}.}
  \item{cexRow}{ See documentation for \link{heatmap}.}
  \item{cexCol}{ See documentation for \link{heatmap}.}
  \item{labRow}{ See documentation for \link{heatmap}.}
  \item{labCol}{ See documentation for \link{heatmap}.}
  \item{main}{ See documentation for \link{heatmap}.}
  \item{xlab}{ See documentation for \link{heatmap}.}
  \item{ylab}{ See documentation for \link{heatmap}.}
  \item{keep.dendro}{ See documentation for \link{heatmap}.}
  \item{verbose}{ See documentation for \link{heatmap}.}
  \item{hExp}{height expansion factor (default is \code{1})}
  \item{wExp}{width expansion factor (default is \code{1})}
  \item{vlines}{vector of positions at which to draw vertical lines} 
  \item{hlines}{vector of positions at which to draw horizontal lines} 
  \item{lasCol}{axis label style (las) for columns}
  \item{\dots}{additional arguments passed along to the image command}
}
\details{
  The new arguments \code{hExp} and \code{wExp} are "expansion factors"
  for the height and width of the figure, respectively. They are used
  to alter the arguments passed to the layout function internally to
  heatmap.

  The new arguments \code{hlines} and \code{vlines} are used to specify
  points at which horizontal or vertical lines, respectively, should be
  overlaid on the image.
}
\value{
  The same as the \link{heatmap} function.
}
\author{
  Kevin R. Coombes \email{krc@silicovore.com}
}
\seealso{
  \code{\link{heatmap}}
}
\examples{
nC <- 30
nR <- 1000
fakeData <- matrix(rnorm(nC*nR), ncol=nC, nrow=nR)
aspectHeatmap(fakeData, scale='none', hExp=2, wExp=1.4, margins=c(6,3))
aspectHeatmap(fakeData, scale='none', hExp=2, wExp=1.4, margins=c(6,3),
              vlines=15.5, hlines=c(100, 300))
}
\keyword{hplot}

