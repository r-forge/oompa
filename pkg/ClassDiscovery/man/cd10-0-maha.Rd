\name{mahalanobisQC}
\alias{mahalanobisQC}
\title{Using Mahalanobis Distance and PCA for Quality Control}
\description{
  Compute the Mahalanobis distance of each sample from the center of an
  \var{N}-dimensional principal component space.
}
\usage{
mahalanobisQC(spca, N)
}
\arguments{
  \item{spca}{object of class \code{\link{SamplePCA}} representing the
    results of a principal components analysis.}
  \item{N}{integer scalar specifying the number of components to use when
    assessing QC.}
}
\details{
  The theory says that, under the null hypothesis that all samples arise
  from the same multivariate normal distribution, the distance from the
  center of a \var{D}-dimensional principal component space should follow a 
  chi-squared distribution with \var{D} degrees of freedom. This theory lets
  us compute p-values associated with the Mahalanobis distances for
  each sample. This method can be used for quality control or outlier
  identification.
}
\value{
  Returns a data frame containing two columns, with the rows
  corresponding to the columns of the original data set on which PCA was
  performed. First column is the chi-squared statistic, with \code{N}
  degrees of freedom. Second column is the associated p-value.
}
\references{
  Coombes KR, et al.\cr
  \emph{Quality control and peak finding for proteomics data collected from
  nipple aspirate fluid by surface-enhanced laser desorption and ionization.}
  Clin Chem 2003; 49:1615-23.
}
\author{
  Kevin R. Coombes \email{krc@silicovore.com}
}
\examples{
library(oompaData)
data(lungData)
spca <- SamplePCA(na.omit(lung.dataset))
mc <- mahalanobisQC(spca, 2)
mc[mc$p.value < 0.01,]
}
\keyword{models}
\keyword{multivariate}

