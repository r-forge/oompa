\name{generics}
\alias{process}
\alias{analyze}
\title{ OOMPA Generic Functions }
\description{
  New generic functions for processing and analyzing microarrays.
}
\usage{
process(object, action, parameter=NULL)
analyze(object, \dots)
}
\arguments{
  \item{object}{ Any OOMPA class representing a microarrays or a set of
    microarrays } 
  \item{action}{ The action to process the class. }
  \item{parameter}{ Any parameters needed to execute the process.}
  \item{\dots}{ Place holder for additional parameters needed in derived
    classes }
}
\details{
  In general, the \code{analyze} method represents an expensive
  computational step carried out in preparation for a graphical
  display, but the semantics may differ from class to class. The default
  implementation of the method performs the null analysis; that is, the
  return value is identical to the object that is passed in as the first
  argument.

  The \code{process} method represents a function that acts on the data
  of some object to process it in some way. For example, normalizing a
  set of microarray data is typically one processing step in a long
  series that is required to take the raw data and turn it into
  something useful.
}
\value{
  Varies depending on the implementation in derived classes.  Typically
  another object of the same or a closely related class.
}
\author{
  Kevin R. Coombes \email{kcoombes@mdanderson.org}
}
\seealso{
  \code{\link{Processor}},
  \code{\link{Pipeline}}
}
\keyword{methods}

