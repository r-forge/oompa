\name{SingleGroup-class}
\docType{class}
\alias{SingleGroup}
\alias{SingleGroup-class}
\alias{as.data.frame,SingleGroup-method}
\alias{plot,SingleGroup,missing-method}
\alias{print,SingleGroup-method}
\alias{summary,SingleGroup-method}

\title{The SingleGroup Class }
\description{  Preliminary analysis of one group of samples for use in
  the \code{\link{SmoothTtest}} class. A key feature is the standard
  quality control plot.}
\usage{
SingleGroup(avg, sd, span = 0.5, name = '')
\S4method{as.data.frame}{SingleGroup}(x, row.names=NULL, optional=FALSE)
\S4method{summary}{SingleGroup}(object, ...)
\S4method{print}{SingleGroup}(x, ...)
\S4method{plot}{SingleGroup,missing}(x, multiple=3, ccl=0, main=x@name,
 xlab='Mean', ylab='Std Dev', xlim=0, ylim=0, ...)
}
\arguments{
  \item{avg}{A numeric vector of mean values}
  \item{sd}{A numeric vector of standard deviations}
  \item{span}{The \code{span} parameter is passed onto \code{\link{loess}}.}
  \item{name}{A character string; the name of this object}
  \item{object}{A \code{SingleGroup} object}
  \item{x}{A \code{SingleGroup} object}
  \item{multiple}{A real number; the multiple of the smoothed standard
    deviation to call significant.}
  \item{ccl}{A list containing objects of the
    \code{\link[PreProcess]{ColorCoding}} class. If left at its default
    value of zero, colors are chosen automatically.}
  \item{main}{Plot title}
  \item{xlab}{Label for the x axis}
  \item{ylab}{Label for the y axis}
  \item{xlim}{Plotting limits for the x axis. If left at the default
    value of zero, then the limits are automatically generated}
  \item{ylim}{Plotting limits for the y axis. If left at the default
    value of zero, then the limits are automatically generated}
  \item{row.names}{See the base version of \code{\link{as.data.frame.default}}}
  \item{optional}{See the base version of \code{\link{as.data.frame.default}}}
  \item{\dots}{The usual extra parameters to generic or plotting routines}
}
\details{
  In 2001 and 2002, Baggerly and Coombes developed the smooth t-test for
  finding differentially expressed genes in microarray data.  Along with
  many others, they began by log-transforming the data as a reasonable
  step in the direction of variance stabilization.  They observed,
  however, that the gene-by-gene standard deviations still seemed to
  vary in a systematic way as a function of the mean log intensity.  By
  borrowing strenght across genes and using \code{\link{loess}} to fit
  the observed standard deviations as a function of the mean, one
  presumably got a better estimate of the true standard deviation.
}
\section{Creating Objects}{
  Objects can be created by calls to the \code{SingleGroup} constructor.
  Users rarely have need to create these objects directly; they are
  usually created as a consequence of the construction of an object of
  the \code{\link{SmoothTtest}} class.
}
\section{Slots}{
  \describe{
    \item{\code{name}:}{A character string; the name of this object}
    \item{\code{avg}:}{The numeric vector of mean values}
    \item{\code{sd}:}{The numeric vector of standard deviations}
    \item{\code{span}:}{The \code{span} parameter used in the
      \code{\link{loess}} function to fit \code{sd} as a function of
      \code{avg}.}
    \item{\code{fit}:}{A list containing components \code{x} and
      \code{y} resulting from the \code{loess} fit.}
    \item{\code{score}:}{A numeric vector; the ratio of the pointwise
      standard deviations to their smooth (loess) estimates.}
  }
}
\section{Methods}{
  \describe{
    \item{as.data.frame(x, row.names=NULL, optional=FALSE)}{Combine the
      slots containing numeric vectors into a data frame, suitable for
      printing or exporting.}
    \item{summary(object, ...)}{Write out a summary of the object.}
    \item{print(x, ...)}{Print the entire object. You never want to do
      this.}
    \item{plot(x, multiple=3, ccl=0, main=x@name, xlab='Mean', ylab='Std
    Dev', xlim=0, ylim=0, ...)}{Produce a scatter plot of the standard
    deviations (\code{x@sd}) as a function of the means (\code{x@avg}).
    The appropriate mutliple of the \code{loess} fit is overlaid, and
    points that exceed this multiple are flaged in a different
    color. Colors in the plotare controlled by the current values of
    \code{\link[PreProcess:semanticColors]{COLOR.CENTRAL.LINE}},
    \code{\link[PreProcess:semanticColors]{COLOR.CONFIDENCE.CURVE}},
    \code{\link[PreProcess:semanticColors]{COLOR.BORING}},
    \code{\link[PreProcess:semanticColors]{COLOR.BAD.REPLICATE}}, and
    \code{\link[PreProcess:semanticColors]{COLOR.WORST.REPLICATE}}.
  }
  }
}
\references{
  Baggerly, K.A., Coombes, K.R., Hess, K.R., Stivers, D.N., Abruzzo,
  L.V., Zhang, W.  Identifying differentially expressed genes in cDNA
  microarray experiments.  J Comp Biol. 8:639-659, 2001.
  
  Coombes, K.R., Highsmith, W.E., Krogmann, T.A., Baggerly, K.A.,
  Stivers, D.N., Abruzzo, L.V.   Identifying and quantifying sources of
  variation in microarray data using high-density cDNA membrane arrays.
  J Comp Biol. 9:655-669, 2002.
}
\author{Kevin R. Coombes <kcoombes@mdanderson.org>}
\seealso{
  \code{\link{SmoothTtest}}
}
\examples{
m <- rnorm(1000, 8, 2.5)
v <- rnorm(1000, 0.7)
plot(m, v)

x <- SingleGroup(m, v, name='bogus')

summary(x)

plot(x)
plot(x, multiple=2)

# cleanup
rm(m, v, x)
}
\keyword{classes}
\keyword{loess}
