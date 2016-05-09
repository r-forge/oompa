\name{palettes}
\alias{alphabet.colors}
\alias{kelly.colors}
\alias{green.armytage.colors}
\alias{palette36.colors}

\title{Polychrome Color Palettes}
\description{
  Four color palettes each containing at least 22 different,
  distinguishable colors.
}
\usage{
alphabet.colors(n=26)
kelly.colors(n=22)
green.armytage.colors(n=26)
palette36.colors(n=36)
}
\arguments{
  \item{n}{an integer; the number of colors desired.}
}
\value{
  Each function returns a character vector of hexadecimal color values
  (such as "#EA9399"). Each color is assigned a name (such as
  "Strong_Pink").  The default value is the maximum number of colors
  available from the individual palette.
}
\author{Kevin R. Coombes <krc@silicovore.com>}
\examples{
palette36.colors(5)
kelly.colors(5)
alphabet.colors(7)
green.armytage.colors(3)
}
\keyword{color}
