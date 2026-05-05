\name{Neighborhood Preservation}
\alias{AverageJaccardDistance}
\alias{AvgContinuity}
\alias{AvgTrustworthiness}
\alias{QNX}

\title{ Metrics for Preservation of Pairwise Distances}
\description{
  Metrics to estimate how well or poorly a dimension reduction algorithm
  preserves local neighborhoods.
}
\usage{
AverageJaccardDistance(M1, M2, k = 50)
AvgContinuity(M1, M2, k = 300)
AvgTrustworthiness(M1, M2, k = 300)
QNX(D1, D2)
}
\arguments{
  \item{M1}{ a data matrix whose rows are the pointsin the original space.}
  \item{M2}{ a data matrix whose rows are the same points after dimension reduction.}
  \item{D1}{ a \code{dist} object containing the distances between
    points in the orignial space.}
  \item{D2}{ a \code{dist} object containing the distances between
    points after dimension reduction.}
  \item{k}{ an integer; the number of neighbors to consider.}
}
\value{
  The \code{Average Jaccard Distance} returns a real number representing the
  score comparing the spaces before and after dimension reduction. The
  otgher three functions return numeric vectors of length \code{k}.
}

\author{Kevin R. Coombes <krc@silicovore.com>}
\examples{
data(prot5) # original space
data(drs) # spaces after diomwnion reduction
AverageJaccardDistance(prot5, drs[["UMAP"]], k = 50)
}
\keyword{htest}
\keyword{models}
