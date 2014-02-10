\name{ColorCoding-class}
\alias{ColorCoding-class}
\alias{ColorCoding}
\alias{colorCode}
\docType{class}
\title{Class "ColorCoding"}
\description{
  A class for associating plotting symbols and colors with a logical
  vector or with levels of a factor.
}
\usage{
ColorCoding(v, color, mark = 1)
colorCode(fac, colorScheme = 1:length(levels(fac)), mark = 1)
}
\arguments{
  \item{v}{ a logical vector }
  \item{color}{ a string or integer representing a color }
  \item{mark}{ an integer representing a plotting symbol, or list of
    plotting symbols that should be associated with levels of the factor
    \code{fac} }
  \item{fac}{ A factor }
  \item{colorScheme}{ A list of colors that should be associated with
    levels of the factor \code{fac}}
}
\details{
  It is often necessary with microarray data to produce multiple plots,
  where each point on the plot corresponds to a gene or a spot on the
  microarray.  Across multiple plots, we often want to use symbols or
  colors to mark subsets of he genes with certain properties.  The
  \code{ColorCoding} class works in tandem with the
  \code{\link{ColorCodedPair}} class to make it easier to maintain
  consistent plotting conventions across multiple graphs.
}
\section{Objects from the Class}{
  Although objects can be created using \code{new}, the preferred method
  is to use the constructor function, \code{ColorCoding}.  To create a
  set of colors and symbols associated with different levels of a
  factor, use the \code{\link{colorCode}} function.
}
\section{Slots}{
  \describe{
    \item{\code{v}}{ a logical vector }
    \item{\code{color}}{ a string or integer representing a color }
    \item{\code{mark}}{ an integer representing a plotting symbol }
  }
}
\section{Methods}{
  There are no specialized methods for objects of this class; all
  activities can be performed by directly accessing the documented slots.
}
\value{
  The constructor returns a valid object of the \code{ColorCoding}
  class.

  The \code{colorCode} function provides a simple interface to associate
  colors and symbols with the levels of a factor. It returns a list of
  valid \code{ColorCoding} objects that can be passed directly to the
  constructor of a \code{\link{ColorCodedPair}} object.
}
\author{
  Kevin R. Coombes \email{krc@silicovore.com},
  P. Roebuck \email{proebuck@mdanderson.org}

}
\seealso{
  \code{\link{ColorCodedPair}}
}
\examples{
showClass("ColorCoding")

theta <- (0:360)*pi/180
x <- cos(theta)
y <- sin(theta)
xp <- x > 0
yp <- y > 0
colors <- list(ColorCoding(xp&yp, oompaColor$BORING),
               ColorCoding(xp&!yp, oompaColor$TOP.TEN),
               ColorCoding(!xp&yp, oompaColor$BOTTOM.TEN),
               ColorCoding(!xp&!yp, oompaColor$CONFIDENCE.CURVE))
plot(ColorCodedPair(x, y, colors))

plot(ColorCodedPair(theta, x, colors))

plot(ColorCodedPair(theta, y, colors),
     xlab='angle in radians', ylab='sine', main='colored sine')

fac <- factor(rep(c('left', 'right'), c(180, 181)))
colors <- colorCode(fac, c('blue', 'red'))
plot(ColorCodedPair(x, y, colors))

## cleanup
rm(x, y, xp, yp, theta, colors, fac)

colorList <- c(oompaColor$BORING, oompaColor$SIGNIFICANT,
               oompaColor$EXPECTED, oompaColor$OBSERVED,
               oompaColor$PERMTEST, oompaColor$FITTED,
               oompaColor$CENTRAL.LINE, oompaColor$CONFIDENCE.CURVE,
               oompaColor$BAD.REPLICATE, oompaColor$WORST.REPLICATE,
               oompaColor$FOLD.DIFFERENCE, oompaColor$BAD.REPLICATE.RATIO,
               oompaColor$TOP.TEN, oompaColor$BOTTOM.TEN,
               oompaColor$TOP.TEN.SOLO, oompaColor$BOTTOM.TEN.SOLO
               )
plot(c(1,4), c(1,4), type='n')
for (i in 1:4) {
  for (j in 1:4) {
    points(i,j, col=colorList[i + 4*(j-1)], pch=16, cex=4)
  }
}
rm(colorList, i, j)
}
\keyword{classes}
\keyword{color}

