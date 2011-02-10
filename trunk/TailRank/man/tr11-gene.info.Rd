\name{gene.info}
\alias{gene.info}
\docType{data}
\title{Gene information for the prostate cancer data set}
\description{
  This data set provides information about the genes included with the
  (partial) prostate cancer data set as part of the tail.rank.test
  package. 
}
\usage{data(gene.info)}
\format{
  A data frame with 2000 observations on the following 6 variables.
  \describe{
    \item{ArrayI.Spot}{a numeric vector; where is this clone spotted on
      the old arrays}
    \item{ArrayII.Spot}{a numeric vector; where is this clone spotted on
    the new arrays}
    \item{Clone.ID}{a factor; the IMAGE clone identifier} 
    \item{Gene.Symbol}{a factor; the official gene symbol}
    \item{Cluster.ID}{a factor; the UniGene cluster number}
    \item{Accession}{a factor; the GenBanlk accession number} 
  }
}
\source{
  The data was originally described in the paper by Lapointe et al., and
  downloaded from the Stanford Microarray Database
  \url{http://genome-www5.stanford.edu/}. We randomly selected 2000 of
  the 42129 spots to include as demonstration data here.
}
\references{
  Lapointe J et al. (2004) Gene expression profiling identifies
  clinically relevant subtypes of prostate cancer. 
  \emph{Proc Natl Acad Sci U S A}, 101, 811--816.
}
\seealso{\code{\link{clinical.info}},
  \code{\link{expression.data}},
  \code{\link{TailRankTest}}
}
\examples{
data(gene.info)
summary(gene.info)
}
\keyword{datasets}
