\name{graph.utility}
\alias{ellipse}
\alias{f.qq}
\alias{f.qt}
\title{OOMPA graphical utility functions}
\description{
  Utility functions for graphics.
}
\usage{
ellipse(a, b, x0=0, y0=0, ...)
f.qq(x, main = "", cut = 0, ...)
f.qt(x, df, main = "", cut = 0, ...)
}
\arguments{
  \item{a}{ Half the length of the elliptical axis in the x-direction }
  \item{b}{ Half the length of the elliptical axis in the y-direction }
  \item{x0}{ X-coordinate of the center of the ellipse}
  \item{y0}{ Y-coordinate of the center of the ellipse}
  \item{main}{ A text string}
  \item{cut}{ A real number}
  \item{df}{ An integer; the number of degrees of freedom in the t-test }
  \item{\dots}{ Additional graphical parameters passed on to lower-level
    functions }
  \item{x}{A numeric vector}
}
\details{
  The \code{ellipse} function draws an ellipse on an existing plots.
  The ellipses produced by this function are oriented with their major
  and minor axes parallel to the coordinate axes. The current
  implementation uses \code{\link{points}} internally.

  The function \code{f.qq} is a wrapper that combines \code{qqnorm} and
  \code{qqline} into a single function call.

  The function \code{f.qt} is a wrapper that produces quantile-quantile
  plots comparing the observed vector \code{x} with a T-distribution.
}
\author{Kevin R. Coombes <kcoombes@mdanderson.org>}

\seealso{ See also \code{\link{points}}}
\examples{
x <- rnorm(1000, 1, 2)
y <- rnorm(1000, 1, 2)
plot(x,y)
ellipse(1, 1, col=6, type='l', lwd=2)
ellipse(3, 2, col=6, type='l', lwd=2)
f.qq(x, main='Demo', col='blue')
f.qq(x, cut = 3)
f.qt(x, df = 3)
f.qt(x, df = 40)
}

\keyword{aplot}
