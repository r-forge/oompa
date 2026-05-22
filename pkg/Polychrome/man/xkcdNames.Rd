\name{xkcdNames}
\alias{xkcdNames}

\title{Consensus Names for Colors}
\description{
  Standardizing color names is a challenging task. Elsewhere in this
  package, we provide tools for assigning either standard names from the
  ISCC or common names from the historical "rgb.txt" files on UNIX
  computers. Here, we provide access toth e "consensus" color names
  derived froma n online survey conducted by Randall Munroe of XKCD.
}
\usage{
xkcdNames(colorset)
}
\arguments{
  \item{colorset}{A character vector containing hexadecimal
    representations of colors.}
}
\details{
  Each of the 949 XKCD consensus color names can be represented by the
  centroid of a region of CIE L*u*v* color space, all of whose points
  should be given the same name. This function takes colors represented
  in the common hexadecimal notation, maps them into L*u*v* color space,
  and assigns the name of the nearest XKCD color.
}
\value{
  A character string containing the standard color name nearest (in CIE
  L*u*v* color space) to each input color.
}
\author{Kevin R. Coombes <krc@silicovore.com>}
\references{
  \url{https://blog.xkcd.com/2010/05/03/color-survey-results/}.
}
\seealso{
  \code{\link{xkcd}}, \code{\link[grDevices]{colors}}.
}
\examples{
data(alphabet)
xkcdNames(alphabet)
colorNames(alphabet)
}
\keyword{color}
