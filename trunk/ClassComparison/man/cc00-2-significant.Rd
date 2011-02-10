\name{significant}
\alias{cutoffSignificant}
\alias{selectSignificant}
\alias{countSignificant}

\title{ Generic Functions for Significance }
\description{
  In the world of multiple testing that is inhabited by most microarray
  or protein profiling experiments, analysts frequently perform separate
  statistical tests for each gene or protein in the experiment.
  Determining cutoffs that achieve statistical significance (in a
  meaningful way) is an inherent part of the procedure.  It is then
  common to select the significant items for further processing or for
  preparing reports, or at least to count the number of significant
  items.  These generic functions provide a standard set of tools for
  selecting and counting the significant items, which can be used with
  various statistical tests and various ways to account for multiple
  testing.
}
\usage{
cutoffSignificant(object, ...)
selectSignificant(object, ...)
countSignificant(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Typically, an object that performs multiple statistical
    tests on microarray or proteomics data.}
  \item{\dots}{All generic methods are designed to take optional
    additional arguments for flexibility in creating derived classes.}
}
\value{
  \code{cutoffSignificant} returns appropriate cutoff values that achieve
  specified signficance criteria.

  \code{selectSignificant} returns a logical vector, with true values
  indicating items that satisfy the cutoff making them statistically
  significant.

  \code{countSignificant} returns an integer, representing the number of
  significant items.
}
\author{Kevin R. Coombes <kcoombes@mdanderson.org>}
\seealso{Classes that implement these methods include
  \code{\link{Bum}}, \code{\link{Sam}}, \code{\link{Dudoit}},
  \code{\link{MultiWilcoxonTest}}, and \code{\link{TNoM}}.
}
\examples{
# Since these are generic methods, there are no examples. See the
# classes that implement specific instances.
}
\keyword{htest}

