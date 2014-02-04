\name{ColorCodedPair-class}
\alias{ColorCodedPair-class}
\alias{ColorCodedPair}
\alias{plot,ColorCodedPair,missing-method}
\docType{class}
\title{Class "ColorCodedPair"}
\description{
  This class represents a vector of (x,y) pairs, each of which should be
  plotted in a specific color with a specific symbol.
}
\usage{
ColorCodedPair(x, y, ccl)
\S4method{plot}{ColorCodedPair,missing}(x, y, \dots)
}
\arguments{
  \item{x}{numeric vector (for \code{ColorCodedPair}) or an
    object of class \code{ColorCodedPair} (for \code{plot})}
  \item{y}{numeric vector}
  \item{ccl}{list of \code{\link{ColorCoding}} objects }
  \item{\dots}{additional arguments are as in the underlying generic methods.}
}
\details{
  It is often necessary with microarray data to produce multiple plots,
  where each point on the plot corresponds to a gene or a spot on the
  microarray. Across multiple plots, we often want to use symbols or colors
  to mark subsets of he genes with certain properties. The
  \code{ColorCodedPair} class works in tandem with the
  \code{\link{ColorCoding}} class to make it easier to maintain consistent
  plotting conventions across multiple graphs.
}
\section{Slots}{
  \describe{
    \item{\code{x}}{numeric vector}
    \item{\code{y}}{numeric vector}
    \item{\code{colorCodingList}}{list of \code{\link{ColorCoding}} objects}
  }
}
\section{Methods}{
  \describe{
    \item{plot(object, \dots)}{Plot the \code{ColorCodedPair} object, with
      appropriate colors and symbols (of course).}
  }
}
\value{
  The constructor returns a valid \code{ColorCodedPair} object.

  The \code{plot} method invisibly returns the object being plotted.
}
\author{
  Kevin R. Coombes \email{krc@silicovore.com},
  P. Roebuck \email{proebuck@mdanderson.org}
}
\seealso{
  \code{\link{colorCode}},
  \code{\link{ColorCoding}}
}
\examples{
showClass("ColorCodedPair")

theta <- (0:360)*pi/180
x <- cos(theta)
y <- sin(theta)
xp <- x > 0
yp <- y > 0
colors <- list(ColorCoding(xp&yp, oompaColor$EXPECTED),
               ColorCoding(xp&!yp, oompaColor$OBSERVED),
               ColorCoding(!xp&yp, oompaColor$PERMTEST),
               ColorCoding(!xp&!yp, oompaColor$FITTED))
plot(ColorCodedPair(x, y, colors))

plot(ColorCodedPair(theta, x, colors))

plot(ColorCodedPair(theta, y, colors),
     xlab='angle in radians', ylab='sine', main='colored sine')

}
\keyword{classes}
\keyword{color}

