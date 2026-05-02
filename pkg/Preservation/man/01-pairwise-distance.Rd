\name{Pairwise Distance Preservation}
\alias{MilnorDistortion}
\alias{SigmaDistortion}
\alias{Stress}
\alias{M1Distortion}
\alias{SpearmanRho}
\alias{EarthMover}

\title{ Metrics for Preservation of Pairwise Distances}
\description{
  Metrics to estimate how well or poorly a dimension reduction algorithm
  preserves pairwise distances between points.
}
\usage{
MilnorDistortion(D1, D2)
SigmaDistortion(D1, D2)
Stress(D1, D2)
M1Distortion(M1, M2)
SpearmanRho(D1, D2)
EarthMover(D1, D2)
}
\arguments{
  \item{M1}{ a data matrix whose rows are the pointsin the original space.}
  \item{M2}{ a data matrix whose rows are the same points after dimension reduction.}
  \item{D1}{ a \code{dist} object containing the distances between
    points in the orignial space.}
  \item{D2}{ a \code{dist} object containing the distances between
    points after dimension reduction.}
}
\value{
  Each of the functions returns a real number representing the
  correpsonding scores comparing the spaces before and after dimension reduction.
}

\author{Kevin R. Coombes <krc@silicovore.com>}
\examples{
data(prot5) # original space
data(drs) # spacers after diomwnion reduction
M1Distortion(prot5, drs[["UMAP"]])
}
\keyword{htest}
\keyword{models}
