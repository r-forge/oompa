\name{colorschemes}
\alias{redscale}
\alias{greenscale}
\alias{bluescale}
\alias{blueyellow}
\alias{cyanyellow}
\alias{redgreen}
\alias{jetColors}
\alias{grayscale}
\alias{greyscale}

\title{Color Schemes for Images and Heat Maps}
\description{
  Create a vector of \code{N} contiguous colors.
}
\usage{
redscale(N)
greenscale(N)
bluescale(N)
blueyellow(N)
cyanyellow(N)
redgreen(N)
jetColors(N)
grayscale(N)
greyscale(N)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{ an integer; the number of distinct levels in the color map }
}
\details{
  The color maps that ship with R (see, for example,
  \code{\link{terrain.colors}}) do not include the most common color
  maps used in publications in the microarray literature. This
  collection of color maps expands the available options.  The functions
  \code{redscale}, \code{greenscale}, and \code{bluescale} each range
  from pure black for low values to a pure primary color for high
  values. The synonyms \code{graysale} and \code{greysale} range from
  pure black to pure white.

  The \code{redgreen} color map ranges from pure green at the low end,
  through black in the middle, to pure red at the high end. Although
  this is the most common color map used in the microarray literature,
  it will prove problematic for individuals with red-green color-blindness. 

  The \code{blueyellow} color map ranges from pure blue at the low end,
  through gray in the middle, to pure yellow at the high end.

  The \code{jetColors} map tries to reproduce the default "jet" color map
  from MATLAB.

  The \code{cyanyellow} color map was added to provide a divergent map
  that should be usable by the majority of individuals whose vision has
  a color deficit. It ranges from cyan (a mixture of blue and green) at
  the low to yellow (a mixture of red and green) at the high end.  Since
  the vast majority of color deficits arise because an individual lacks
  cones for one of the three primary colors (red, green, or blue), this
  colormap should stil provide adequate contrasts.
}
\value{
  A character vector '\code{cv}' of color names. This can be used to
  create a user-defined color palette for subsequent graphics by
  '\code{palette(cv)}' or directly in a '\code{col=}' specification in
  '\code{par}' or in graphics functions such as '\code{image}' or
  '\code{heatmap}'. 
}
\section{BUGS}{
  The names \code{redgreen} and \code{blueyellow} are inconsistent with
  respect to which color represents low values and which color
  represents high values. It is too late to fix this.
}
\author{Kevin R. Coombes <krc@silicovore.com>}
\seealso{\code{\link{rainbow}}, \code{\link{topo.colors}},
  \code{\link{terrain.colors}}, \code{\link{heat.colors}},
  \code{\link{rgb}}, \code{\link[graphics]{image}}, \code{\link{heatmap}},
  \code{\link{palette}}.
}
\examples{
data <- matrix(1:1024, nrow=1024)
image(data, col=bluescale(64))
image(data, col=redgreen(32))
image(data, col=redscale(128))
image(data, col=blueyellow(64))
image(data, col=jetColors(64))
image(data, col=grayscale(64))
image(data, col=cyanyellow(64))
rm(data) # cleanup
}
\keyword{ color }

