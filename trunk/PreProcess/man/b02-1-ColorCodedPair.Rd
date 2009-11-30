\name{ColorCodedPair-class}
\docType{class}
\alias{ColorCodedPair-class}
\alias{ColorCodedPair}
\alias{plot,ColorCodedPair,missing-method}
\title{ The ColorCodedPair class }
\description{
  This class represents a vector of (x,y) pairs, each of which should be
  plotted in a specific color with a specific symbol.
}
\usage{
ColorCodedPair(x, y, ccl)
\S4method{plot}{ColorCodedPair,missing}(x, y, ...)
}
\arguments{
  \item{x}{ A numeric vector (for \code{ColorCodedPair}) or a
    ColorCodedPair object (for \code{plot})}
  \item{y}{ A numeric vector }
  \item{ccl}{ A list of \code{\link{ColorCoding}} objects }
  \item{\dots}{Additional arguments are as in the underlying generic methods.}
}
\details{
  It is often necessary with microarray data to produce multiple plots,
  where each point on the plot corresponds to a gene or a spot on the
  microarray.  Across multiple plots, we often want to use symbols or
  colors to mark subsets of he genes with certain properties.  The
  \code{ColorCodedPair} class works in tandem with the
  \code{\link{ColorCoding}} class to make it easier to maintain
  consistent plotting conventions across multiple graphs.
}
\section{Slots}{
  \describe{
    \item{\code{x}}{ A numeric vector }
    \item{\code{y}}{ A numeric vector }
    \item{\code{colorCodingList}}{ A list of \code{\link{ColorCoding}}
      objects }
  }
}
\section{Methods}{
  \describe{
    \item{plot(object, ...)}{Plot the \code{ColorCodedPair} object, with
      appropriate colors and symbols (of course).}
  }
}
\value{
  The constructor returns a valid \code{ColorCodedPair} object.

  The \code{plot} method invisibly returns the object being plotted.
}
\author{Kevin R. Coombes <kcoombes@mdanderson.org>}
\seealso{The \code{\link{ColorCoding}} class, \code{\link{colorCode}}}
\examples{
theta <- (0:360)*pi/180
x <- cos(theta)
y <- sin(theta)
xp <- x > 0
yp <- y > 0
colors <- list(ColorCoding(xp&yp, COLOR.EXPECTED),
               ColorCoding(xp&!yp, COLOR.OBSERVED),
               ColorCoding(!xp&yp, COLOR.PERMTEST),
               ColorCoding(!xp&!yp, COLOR.FITTED))
plot(ColorCodedPair(x, y, colors))

plot(ColorCodedPair(theta, x, colors))

plot(ColorCodedPair(theta, y, colors),
     xlab='angle in radians', ylab='sine', main='colored sine')

}
\keyword{color}

