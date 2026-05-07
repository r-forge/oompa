\name{More Path Preservation Metrics}
\alias{Smoothness}
\alias{RankSmoothness}
\alias{Coil}
\alias{Contact}
\alias{EndpointDisplacement}
\alias{CrossingNumber}
\title{ More Metrics for Path Preservation}
\description{
  Metrics to estimate how well or poorly a dimension reduction algorithm
  preserves paths
}
\usage{
Smoothness(M1, M2)
RankSmoothness(M1, M2)
Coil(M1, M2)
Contact(M1, M2, q = 0.1)
EndpointDisplacement(M1, M2, k = 5)
CrossingNumber(M2, eps = 1e-8)
}
\arguments{
  \item{M1}{ a data matrix whose rows are the pointsin the original space.}
  \item{M2}{ a data matrix whose rows are the same points after
    dimension reduction.}
  \item{q}{a real number; the quantile of pairwise diatances used to define "close".}
  \item{k}{an intefger; the number of points to average at each end of
    the path.}
  \item{eps}{A real fudge factor (defauly 1E-8) used when normalizing
    the crossing number.}
}
\value{
  Each of the functions returns a real number representing the
  corresponding scores comparing the spaces before and after dimension reduction.
}

\author{Kevin R. Coombes <krc@silicovore.com>}
\examples{
data(prot5) # original space
data(drs) # spaces after diomwnion reduction
Smoothness(prot5, drs[["TSNE"]])
}
\keyword{htest}
\keyword{models}
