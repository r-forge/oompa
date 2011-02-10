\name{clinical.info}
\alias{clinical.info}
\docType{data}
\title{Experimental info for the prostate cancer data set}
\description{
  This data set provides experimental and clinical information about the
  (partial) prostate cancer data set included for demonstration purposes
  as part of the tail.rank.test package. The experiments were two-color
  glass microarrays printed at Stanford.
}
\usage{data(clinical.info)}
\format{
  A data frame with 112 observations on the following 6 variables.
  \describe{
    \item{Arrays}{A factor containing the barcode of the microarray on
      which the experiment was performed. Each of the 112 entries should
    be distinct.}
    \item{Reference}{A factor describing the reference sample used in
      each experiment. This was a common reference, so the identifiers
      here are not meaningful.}
    \item{Sample}{A factor identifying the test sample in each
      experiment. These match the codes published in the original
      paper.}
    \item{Status}{A factor with three levels identifying normal prostate
      (\code{N}), prostate cancer (\code{T}), or lymph node metastasis
      (\code{L}).}
    \item{Subgroups}{A factor with five levels: \code{I} \code{II}
      \code{III} \code{N} \code{O}. These correspond to the groups found
      in the original paper using clustering.}
    \item{ChipType}{a factor with levels \code{new} or \code{old}. At
      least two different print designs of microarrays were used in this
    experiment; this factor identifies the design.}
  }
}
\source{
  The data was originally described in the paper by Lapointe et al., and
  downloaded from the Stanford Microarray Database
  \url{http://genome-www5.stanford.edu/}.
}
\references{
  Lapointe J et al. (2004) Gene expression profiling identifies
  clinically relevant subtypes of prostate cancer. 
  \emph{Proc Natl Acad Sci U S A}, 101, 811--816.
}
\seealso{\code{\link{expression.data}},
  \code{\link{gene.info}},
  \code{\link{TailRankTest}}
}
\examples{
data(clinical.info)
summary(clinical.info) 
}
\keyword{datasets}
