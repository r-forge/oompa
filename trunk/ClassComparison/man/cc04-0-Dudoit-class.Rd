\name{Dudoit-class}
\docType{class}
\alias{Dudoit}
\alias{Dudoit-class}
\alias{plot,Dudoit,missing-method}
\alias{cutoffSignificant,Dudoit-method}
\alias{selectSignificant,Dudoit-method}
\alias{countSignificant,Dudoit-method}

\title{The Dudoit Class }
\description{
  An implementation of the method of Dudoit and colleagues to apply the
  Westfall-Young adjustment to p-values to control the family-wise error
  rate when analyzing microarray data.
}
\usage{
Dudoit(data, classes, nPerm = 1000, verbose = TRUE)
\S4method{plot}{Dudoit,missing}(x, y, xlab='T-Statistic', ylab='P-Value', ...)
\S4method{cutoffSignificant}{Dudoit}(object, alpha, ...)
\S4method{selectSignificant}{Dudoit}(object, alpha, ...)
\S4method{countSignificant}{Dudoit}(object, alpha, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ Either a data frame or matrix with numeric values or an
    \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} as defined
    in the BioConductor tools for analyzing microarray data.} 
  \item{classes}{ If \code{data} is a data frame or matrix, then classes
    must be either a logical vector or a factor. If \code{data} is an
    \code{ExpressionSet}, then \code{classes} can be a character string that
    names one of the factor columns in the associated
    \code{\link[Biobase]{phenoData}} subobject.}
  \item{nPerm}{An integer; the number of permutations to perform}
  \item{verbose}{A logical flag}
  \item{object}{A \code{Dudoit} object}
  \item{alpha}{A real number; the target family-wise error rate}
  \item{x}{A \code{Dudoit} object}
  \item{y}{Nothing, since it is supposed to be mssing. Changes to the Rd
    processor require docuemnting the missing entry.}
  \item{xlab}{Label for the x axis}
  \item{ylab}{Label for the y axis}
  \item{\dots}{The usual extra arguments for generic or plotting routines.}
}
\details{
  In 2002, Dudoit and colleagues introduced a method to adjust the
  p-values when performing gene-by-gene tests for differential
  expression.  The adjustment was based on the method of Westfall and
  Young, with the goal of controlling the family-wise error rate.
}
\section{Objects from the Class}{
  As usual, objects can be created by \code{new}, but better methods are
  available in the form of the \code{Dudoit} function. The basic
  inputs to this function are the same as those used for row-by-row
  statistical tests throughout the ClassComparison package; a detailed
  description can be found in the \code{\link{MultiTtest}} class.

  The additional input determines the number, \code{nPerm}, of
  permutations to perform. The accuracy of the p-value adjustment
  depends on this value.  Since the implementation is in R (and does not
  call out to something compiled like C or FORTRAN), however, the
  computations are slow. The default value of 1000 can take a long
  time with modern microarrays that contain 40,000 spots.
}
\section{Slots}{
  \describe{
    \item{\code{adjusted.p}:}{The numeric vector of adjusted p-values.}
    \item{\code{t.statistics}:}{Object of class \code{numeric}
      containing the computed t-statistics. }
    \item{\code{p.values}:}{Object of class \code{numeric} containing
      the computed p-values. }
    \item{\code{groups}:}{Object of class \code{character} containing
      the names of the classes being compared. }
    \item{\code{call}:}{Object of class \code{call} containing the
      function call that created the object. }
  }
}
\section{Extends}{
Class \code{MultiTtest}, directly.  In particular, objects of this class
inherit methods for \code{summary}, \code{hist}, and \code{plot} from
the base class.
}
\section{Methods}{
  \describe{
    \item{cutoffSignificant(object, alpha, ...)}{Determine cutoffs on
      the adjusted p-values at the desired significance level.  In other
      words, this function simplyt returns \code{alpha}.}
    \item{selectSignificant(object, alpha, ...)}{Compute a logical
      vector for selecting significant test results.}
    \item{countSignificant(object, alpha, ...)}{Count the number of
      significant test results.}
    \item{plot}{\code{signature(x = Dudoit, y = missing)}: ... }
  }
}
\value{
  The standard method for \code{plot} returns what you would expect.

  The \code{cutoffSignificant} method returns a real number (its input
  value \code{alpha}).  The \code{selectSignficant} method returns a
  vector of logical values identifying the significant test results, and
  \code{countSignificant} returns an integer counting the number of
  significant test results.
}
\references{
  Dudoit, S., Y.H. Yang, M.J. Callow, and T.P. Speed. 2002. Statistical
  Methods for Identifying Differentially Expressed Genes in Replicated
  cDNA Microarray Experiments, Statistica Sinica, 12(1): 111-139.

  Westfall, P.H., Young, S.S. Resampling-based multiple testing: examples
  and methods for p-value adjustment.  Wiley series in probability and
  mathematics statistics. John Wiley and Sons, 1993.
}
\author{Kevin R. Coombes <kcoombes@mdanderson.org>}
\seealso{
  \code{\link{MultiTtest}}, \code{\link{Bum}}, \code{\link{SmoothTtest}}
}
\examples{
ng <- 10000
ns <- 15
nd <- 200
fake.class <- factor(rep(c('A', 'B'), each=ns))
fake.data <- matrix(rnorm(ng*ns*2), nrow=ng, ncol=2*ns)
fake.data[1:nd, 1:ns] <- fake.data[1:nd, 1:ns] + 2
fake.data[(nd+1):(2*nd), 1:ns] <- fake.data[(nd+1):(2*nd), 1:ns] - 2

# the permutation test is slow. it really needs more than
# 100 permutations, but this is just an example...
dud <- Dudoit(fake.data, fake.class, nPerm=100)
summary(dud)
plot(dud)
countSignificant(dud, 0.05)

rm(ng, ns, nd, fake.class, fake.data, dud)
}
\keyword{classes}
\keyword{multivariate}
\keyword{htest}
