\name{path Preservation}
\alias{LengthDistortion}
\alias{SegmentVariance}
\alias{Curvature}
\alias{SpatialSimilarity}

\title{ Metrics for Path Preservation }
\description{
  Utility functions for computing vectors of row-by-row means,
  variances, and t-statistics.
}
\usage{
LengthDistortion(M1, M2)
SegmentVariance(M1, M2)
Curvature(M1, M2)
SpatialSimilarity(M1, M2)
}
\arguments{
  \item{M1}{ a data matrix decribing the path in the original space.}
  \item{M2}{ a data matrix describing the (same) path after dimension reduction.}
}
\value{
  Each of the four functions returns a real number representing the log
  ratio of the scores betweent he two path representation.
}

\author{Kevin R. Coombes <krc@silicovore.com>}
\examples{
}
\keyword{htest}
