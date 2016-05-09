\name{isccNames}
\alias{isccNames}

\title{Standard Names for Colors}
\description{
  The Inter-Society Color Council, in cooperation with the United States
  National Bureau of Standards, developed a list of 267 standardized
  color names.  We provide tools to convert hexadecimal colors to
  standardized names. 
}
\usage{
isccNames(colorset)
}
\arguments{
  \item{colorset}{A character vector containing hexadecimal
    representations of colors.}
}
\details{
  Each of the ISCC-NBS 267 standard color names is represented by the
  centroid of a region of CIE L*u*v* color space, all of whose points
  should be given the same name.  Colors represented in the common
  hexadecimal notation are mapped into L*u*v* color space and assigned
  the name of the nearest centroid. 
}
\value{
  A character string containing the standard color name nearest (in CIE
  L*u*v* color space) to each input color.
}
\author{Kevin R. Coombes <krc@silicovore.com>}
\references{
  Kelly KL. Twenty-Two Colors of Maximum Contrast. Color Eng., 1965;
  3:26--7. 

  Also see the Inter-Society Color Council web site
  (\url{http://www.iscc.org/}) and the Texas
  Precancel Club (\url{http://tx4.us/nbs-iscc.htm}). 
}
\seealso{
  \code{\link{iscc}}
}
\examples{
data(alphabet)
isccNames(alphabet)
}
\keyword{color}
