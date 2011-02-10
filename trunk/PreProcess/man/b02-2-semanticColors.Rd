\name{semantic colors}
\alias{semanticColors}
\alias{COLOR.BORING}
\alias{COLOR.SIGNIFICANT}
\alias{COLOR.EXPECTED}
\alias{COLOR.OBSERVED}
\alias{COLOR.PERMTEST}
\alias{COLOR.FITTED}
\alias{COLOR.CENTRAL.LINE}
\alias{COLOR.CONFIDENCE.CURVE}
\alias{COLOR.BAD.REPLICATE}
\alias{COLOR.WORST.REPLICATE}
\alias{COLOR.FOLD.DIFFERENCE}
\alias{COLOR.BAD.REPLICATE.RATIO}
\alias{COLOR.TOP.TEN}
\alias{COLOR.BOTTOM.TEN}
\alias{COLOR.BOTTOM.TEN.SOLO}
\alias{COLOR.TOP.TEN.SOLO}

\docType{data}
\title{ Pre-defined colors for microarray plots }
\description{
  A collection of predefined color names to help ensure consistency in
  multiple graphical displays of microarray data.

  \describe{
    \item{\code{COLOR.BORING}:}{Used to mark uninteresting points in a
      plot; gray.}
    \item{\code{COLOR.SIGNIFICANT}:}{Used to mark points that are
      statistically significant; red}
    \item{\code{COLOR.EXPECTED}:}{Used to draw curves representing
      an expected distribution; blue}
    \item{\code{COLOR.OBSERVED}:}{Used to draw curves indicating the
      observed distribution; darkgreen}
    \item{\code{COLOR.PERMTEST}:}{Used to draw curves indicating
      distributions derived from a permutation test; magenta}
    \item{\code{COLOR.FITTED}:}{Used to draw curves obtained by some
      fitting procedure, such as loess; orange}
    \item{\code{COLOR.CENTRAL.LINE}:}{Used to draw lines through the
      centers of distributions or expected values; blue}
    \item{\code{COLOR.CONFIDENCE.CURVE}:}{Used to draw confidence bounds
      around curves; red3}
    \item{\code{COLOR.BAD.REPLICATE}:}{Used to indicate highly variable
      points; purple1}
    \item{\code{COLOR.WORST.REPLICATE}:}{Used to mark extraordinarily
      variable points; purple3}
    \item{\code{COLOR.FOLD.DIFFERENCE}:}{Used to indicate points with
      large fold difference; skyblue}
    \item{\code{COLOR.BAD.REPLICATE.RATIO}:}{Used to flag points for
      which the ratios at replicate spots are highly variable; violetred}
    \item{\code{COLOR.TOP.TEN}:}{Used to mark points in the "top ten"
      list; cadetblue}
    \item{\code{COLOR.BOTTOM.TEN}:}{Used to mark points in "bottom ten"
      list of most underexpressed genes; pink}
    \item{\code{COLOR.BOTTOM.TEN.SOLO}:}{Use unknown; palegreen}
    \item{\code{COLOR.TOP.TEN.SOLO}:}{Use unknown; deeppink}
  }
}
\examples{
x <- seq(0, 2*pi, by=0.1)
plot(x, sin(x), col=COLOR.BORING)
}
\keyword{datasets}
