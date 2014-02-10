\name{semantic colors}
\alias{oompaColor}
\title{ Pre-defined colors for microarray plots }
\description{
  A collection of predefined color names to help ensure consistency in
  multiple graphical displays of microarray data.

  \describe{
    \item{\code{oompaColor}:}{a list containing named components
      allowing the user to systematically use colors for different
    interpretations.}
    \item{\code{oompaColor$BORING}:}{Used to mark uninteresting points in a
	plot; gray.}
    \item{\code{oompaColor$SIGNIFICANT}:}{Used to mark points that are
      statistically significant; red}
    \item{\code{oompaColor$EXPECTED}:}{Used to draw curves representing
      an expected distribution; blue}
    \item{\code{oompaColor$OBSERVED}:}{Used to draw curves indicating the
      observed distribution; darkgreen}
    \item{\code{oompaColor$PERMTEST}:}{Used to draw curves indicating
      distributions derived from a permutation test; magenta}
    \item{\code{oompaColor$FITTED}:}{Used to draw curves obtained by some
      fitting procedure, such as loess; orange}
    \item{\code{oompaColor$CENTRAL.LINE}:}{Used to draw lines through the
      centers of distributions or expected values; blue}
    \item{\code{oompaColor$CONFIDENCE.CURVE}:}{Used to draw confidence bounds
      around curves; red3}
    \item{\code{oompaColor$BAD.REPLICATE}:}{Used to indicate highly variable
      points; purple1}
    \item{\code{oompaColor$WORST.REPLICATE}:}{Used to mark extraordinarily
      variable points; purple3}
    \item{\code{oompaColor$FOLD.DIFFERENCE}:}{Used to indicate points with
      large fold difference; skyblue}
    \item{\code{oompaColor$BAD.REPLICATE.RATIO}:}{Used to flag points for
      which the ratios at replicate spots are highly variable; violetred}
    \item{\code{oompaColor$TOP.TEN}:}{Used to mark points in the "top ten"
      list; cadetblue}
    \item{\code{oompaColor$BOTTOM.TEN}:}{Used to mark points in "bottom ten"
      list of most underexpressed genes; pink}
    \item{\code{oompaColor$BOTTOM.TEN.SOLO}:}{Use unknown; palegreen}
    \item{\code{oompaColor$TOP.TEN.SOLO}:}{Use unknown; deeppink}
  }
}
\examples{
oompaColor
x <- seq(0, 2*pi, by=0.1)
plot(x, sin(x), col=oompaColor$BORING)
}
\keyword{color}

