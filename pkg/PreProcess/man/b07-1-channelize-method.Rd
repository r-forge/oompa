\name{channelize-method}
\alias{channelize}
\alias{channelize-method}
\alias{channelize,ANY-method}
\docType{methods}
\title{Method "channelize"}
\description{
  \code{channelize} is a generic function used to propagate the class of
  derived objects through a processing pipeline.
}
\usage{
\S4method{channelize}{ANY}(object, \dots)
}
\arguments{
  \item{object}{an object for which pipeline propagation is desired}
  \item{\dots}{additional arguments affecting the elapsed time produced}
}
\details{
  Having abstracted away the notion of extracting a particular
  measurement from a \code{\link{CompleteChannel}} object and producing 
  a simple \code{Channel}, we need a way to allow object-oriented
  programming and derived classes to work with our
  \code{\link{Processor}} and \code{\link{Pipeline}} routines. The
  underlying idea is that specific kinds of microarrays or specific
  software to quantify microarrays might have special properties that
  should be exploited in processing. For example, the first few
  generations of microarrays printed at M.D. Anderson spotted every cDNA
  clone in duplicate. The analysis of such arrays should exploit this
  additional structure. In order to do so, we must derive classes from
  \code{CompleteChannel} and \code{Channel} and ensure that the classes
  of extracted objects are propagated correctly through the processing
  pipeline. The \code{channelize} method achieves this goal.
}
\value{
  Returns a string, which represents the name of a class (suitable for
  passing to the \code{new} constructor) extracted from an object
  belonging to a class derived from \code{\link{CompleteChannel}}.
}
\note{
  The sections above document the method's usage by OOMPA's pipeline,
  not the actual intent of the generic itself.
}
\author{
  Kevin R. Coombes \email{krc@silicovore.com},
  P. Roebuck \email{proebuck@mdanderson.org}
}
\seealso{
  \code{\link{Channel}},
  \code{\link{CompleteChannel}},
  \code{\link{Pipeline}},
  \code{\link{Processor}}
}
\keyword{methods}
\keyword{manip}

