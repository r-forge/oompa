\name{channelize}
\alias{channelize}
\title{ Object-oriented conversion from CompleteChannel to Channel}
\description{
  Provide a generic function for propagating the class of derived
  objects through a processing pipeline.
}
\usage{
channelize(object)
}
\arguments{
  \item{object}{ An object of a clas derived from \code{\link{CompleteChannel}}.}
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
  additional structure.  In order to do so, we must derive classes from
  \code{CompleteChannel} and \code{Channel} and ensure that the classes
  of extracted objects are propagated correctly through the processing
  pipeline. The \code{channelize} methods achieves this goal.

}
\value{
  Returns a string, which represents the name of a class (suitable for
  passing to the \code{new} constructor) extracted from an object
  belonging to a class derived from \code{\link{CompleteChannel}}.
}
\author{
  Kevin R. Coombes \email{kcoombes@mdanderson.org}
}
\seealso{
  \code{\link{Channel}},
  \code{\link{CompleteChannel}},
  \code{\link{Pipeline}},
  \code{\link{Processor}}
}
\keyword{methods}
\keyword{manip}

