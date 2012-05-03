\name{generics}
\alias{process}
\alias{process,ANY-method}
\alias{analyze}
\alias{analyze,ANY-method}
\title{Methods "process" and "analyze"}
\description{
  New generic functions for processing and analyzing microarrays.
}
\usage{
\S4method{process}{ANY}(object, action, parameter=NULL)
\S4method{analyze}{ANY}(object, \dots)
}
\arguments{
  \item{object}{any OOMPA class representing a microarrays or a set of
    microarrays} 
  \item{action}{the action to process the class}
  \item{parameter}{any parameters needed to execute the process}
  \item{\dots}{extra arguments for generic routines}
}
\details{
  In general, the \code{analyze} method represents an expensive computational
  step carried out in preparation for a graphical display, but the semantics
  may differ from class to class. The default implementation of the method
  performs the null analysis; that is, the return value is identical to the
  object that is passed in as the first argument.

  The \code{process} method represents a function that acts on the data
  of some object to process it in some way. For example, normalizing a
  set of microarray data is typically one processing step in a long series
  that is required to take the raw data and turn it into something useful.
}
\value{
  The form of the value returned by either \code{process} or \code{analyze}
  depends on the class of its argument. See the documentation of the particular
  methods for details of what is produced by that method.
}
\author{
  Kevin R. Coombes \email{kcoombes@mdanderson.org},
  P. Roebuck \email{proebuck@mdanderson.org}
}
\seealso{
  \code{\link{Pipeline}},
  \code{\link{Processor}}
}
\keyword{methods}

