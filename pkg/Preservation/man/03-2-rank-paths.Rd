\name{Path RankPreservation}
\alias{RankLengthDistortion}
\alias{RankCurvature}
\alias{RankSpatialSimilarity}

\title{ Rank-Based Metrics for Path Preservation }
\description{
  Rank-based metrics to estimate how well or poorly a dimension
  reduction algorithm preserves pathss.
}
\usage{
RankLengthDistortion(M1, M2)
RankCurvature(M1, M2)
RankSpatialSimilarity(M1, M2)
}
\arguments{
  \item{M1}{ a data matrix decribing the path in the original space.}
  \item{M2}{ a data matrix describing the (same) path after dimension reduction.}
}
\value{
  Each of the functions returns a real number representing the log
  ratio of the scores between the two path representation.
}

\author{Kevin R. Coombes <krc@silicovore.com>}
\examples{
data(prot5) # original space
data(drs) # spaces after diomwnion reduction
RankCurvature(prot5, drs[["DM"]])
}
\keyword{htest}
\keyword{models}
