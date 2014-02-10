\name{expression.data}
\alias{expression.data}
\docType{data}
\title{Microarray expression data on prostate cancer}
\description{
  A subset of the microarray data from a study of prostate
  cancer at Stanford is supplied as demo data with the
  tail.rank.test package.
}
\usage{data(expression.data)}
\format{
  A data frame with 2000 observations on the 112 variables.
  Each column represent a different patient sample, as
  described in the accompanying data.frame called
  \code{\link{clinical.info}}.
}
\details{
  This data set contains normalized microarray expression
  data on 2000 randomly selected genes from a prostate
  cancer data set. The study was originially described
  in a publication by Lapointe et al. The experiments
  were performed on two-color glass microarrays printed
  at Stanford and available from the Stanford Microarray
  Database.  We downloaded the raw data and preprocessed
  it. In particular,after background correction and loess
  normalization, we computed log ratios between the channels.
  We then randomly selected 2000 of the 42129 spots to
  include as demonstration data here.
}
\source{
  \url{http://genome-www5.stanford.edu/}.
}
\references{
  Lapointe J et al. (2004) Gene expression profiling identifies
  clinically relevant subtypes of prostate cancer. 
  \emph{Proc Natl Acad Sci U S A}, 101, 811--816.
}
\seealso{\code{\link{clinical.info}},
  \code{\link{gene.info}}
}
\examples{
data(expression.data)
summary(expression.data)
}
\keyword{datasets}
