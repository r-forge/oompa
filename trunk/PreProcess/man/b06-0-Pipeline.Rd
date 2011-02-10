\name{Pipeline-class}
\docType{class}
\alias{Pipeline}
\alias{Pipeline-class}
\alias{process,ANY,Pipeline-method}
\alias{summary,Pipeline-method}
\alias{PIPELINE.STANDARD}
\alias{PIPELINE.MDACC.DEFAULT}
\alias{makeDefaultPipeline}

\title{ The Pipeline class }
\description{
  A \code{Pipeline} represents a standard multi-step procedure for
  processing microarray data.  A \code{Pipeline} represents a series of
  \code{\link{Processor}}s that should be applied in order. You can
  think of a pipeline as a completely defined (and reusable) set of
  transformations that is applied uniformly to every microarray in a
  data set.
}
\usage{
\S4method{process}{ANY,Pipeline}(object, action, parameter=NULL)
\S4method{summary}{Pipeline}(object, ...)
makeDefaultPipeline(ef = PROC.SIGNAL, ep = 0,
                    nf = PROC.GLOBAL.NORMALIZATION, np = 0,
                    tf = PROC.THRESHOLD, tp = 25,
                    lf = PROC.LOG.TRANSFORM, lp = 2,
                    name = "standard pipe",
                    description = "my method")
}
\arguments{
  \item{object}{In the \code{process} method, any object appropriate for
    the input to the \code{Pipeline}. In the \code{summary} method, a
    \code{Pipeline} object. } 
  \item{action}{ A \code{Pipeline} object used to process an object. }
  \item{parameter}{ Irrelevant, since the \code{Pipline} ignores the
    parameter when \code{process} is invoked.}
  \item{\dots}{Additional arguments are as in the underlying generic methods.}
  \item{ef}{ ``Extractor function'': First \code{\link{Processor}} in
    the \code{Pipeline}, typically a method that extracts a single kind
    of raw measurement from a microarray}
  \item{ep}{ Default parameter value for \code{ef}}
  \item{nf}{ ``Normalization function'' : Second \code{\link{Processor}} in
    the \code{Pipeline}, typically a normalization step. }
  \item{np}{ Default parameter value for \code{nf}}
  \item{tf}{ ``Threshold function'' : Third \code{\link{Processor}} in
    the \code{Pipeline}, typically a step that truncates data below at
    some threshold. }
  \item{tp}{ Default parameter value for \code{tf}}
  \item{lf}{ ``Log function'' : Fourth \code{\link{Processor}} in
    the \code{Pipeline}, typically a log transformation. }
  \item{lp}{ Default parameter value for \code{lf}}
  \item{name}{A string; the name of the pipeline}
  \item{description}{A string; a longer description of the pipeline}
}
\details{
  A key feature of a \code{Pipeline} is that it is supposed to represent
  a standard algorithm that is applied to all objects when processing a
  microarray data set. For that reason, the \code{parameter} that can be
  passed to the \code{process} function is ignored, ensuring that the
  same parameter values are used to process all objects. By contrast,
  each \code{\link{Processor}} that is inserted into a \code{Pipeline}
  allows the user to supply a parameter that overrides its default
  value.

  We provide a single constructor, \code{makeDefaultPipeline} to build a
  specialized kind of \code{Pipeline}, tailored to the analysis of
  fluorescently labeled single channels in a microarray experiment. More
  general \code{Pipeline}s can be constructed using \code{new}.
}
\section{Slots}{
  \describe{
    \item{\code{proclist}:}{A list of \code{\link{Processor}} objects.}
    \item{\code{name}:}{A string containing the name of the object}
    \item{\code{description}:}{A string containing a longer description
      of the object}
  }
}
\section{Methods}{
  \describe{
    \item{process(object, action, parameter)}{Apply the series of
      functions represented by the \code{Pipeline} \code{action} to the
      object, updating its history appropriately. The \code{parameter}
      is ignored, since the \code{Pipeline} always uses its default
      values.
    }
    \item{summary(object, ...)}{Write out a summary of the object.}
  }
}
\value{
  The return value of the generic function \code{process} is always
  an object related to its input, which keeps a record of its
  history. The precise class of the result depends on the functions used
  to create the \code{Pipeline}. 
}
\section{Pre-defined Pipelines}{
  The library comes with two \code{Pipeline} objects already defined
  \describe{
    \item{\code{PIPELINE.STANDARD}}{Takes a \code{\link{Channel}} object
      as input. Performs global normalization by rescaling the 75th
      percentile to 1000, truncates below at 25, then performs log
      (base-two) transformation.
    }
    \item{\code{PIPELINE.MDACC.DEFAULT}}{Takes a
      \code{\link{CompleteChannel}} as input, extracts the raw signal
      intensity, and then performs the same processing as
      \code{PIPELINE.STANDARD}. 
    }
  }
}
\author{Kevin R. Coombes <kcoombes@mdanderson.org>}
\seealso{ \code{\link{Channel}}, \code{\link{process}},
  \code{\link{CompleteChannel}}
}
\examples{
# simulate a moderately realistic looking microarray
nc <- 100
nr <- 100
v <- rexp(nc*nr, 1/1000)
b <- rnorm(nc*nr, 80, 10)
s <- sapply(v-b, max, 1)
ct <- ChannelType('user', 'random', nc, nr,  'fake')
subbed <- Channel(name='fraud', parent='', type=ct, vec=s)
rm(ct, nc, nr, v, b, s)		# clean some stuff

# example of standard data processing
processed <- process(subbed, PIPELINE.STANDARD)

summary(processed)

par(mfrow=c(2,1))
plot(processed)
hist(processed)

par(mfrow=c(1,1))
image(processed)

rm(subbed, processed)
}
\keyword{manip}
