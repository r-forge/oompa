\name{Processor-class}
\docType{class}
\alias{Processor}
\alias{Processor-class}
\alias{summary,Processor-method}
\alias{process,Processor-method}
\alias{process,Channel,Processor-method}
\alias{PROC.SUBTRACTOR}
\alias{PROC.THRESHOLD}
\alias{PROC.GLOBAL.NORMALIZATION}
\alias{PROC.LOG.TRANSFORM}
\alias{PROC.MEDIAN.EXPRESSED.NORMALIZATION}
\alias{PROC.SUBSET.NORMALIZATION}
\alias{PROC.SUBSET.MEAN.NORMALIZATION}

\title{ The Processor class }
\description{
  A \code{Processor} represents a function that acts on the data of a
  some object to process it in some way. The result is always
  another related object, which should record some history about
  exactly how it was processed.
}
\usage{
\S4method{process}{Channel,Processor}(object, action, parameter=NULL)
\S4method{summary}{Processor}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{In the \code{process} method, a \code{Channel}
    object. In the \code{summary} method, a \code{Processor} object}
  \item{action}{ A \code{Processor} object used to process a
    \code{Channel}. }
  \item{parameter}{ Any object that makes sense as a parameter to the
    function represented by the \code{Processor} \code{action}}
  \item{\dots}{Additional arguments are as in the underlying generic methods.}
}
\section{Slots}{
  \describe{
    \item{\code{f}:}{A function that will be used to process
      microarray-related object}
    \item{\code{default}:}{The default value of the parameters to the
      function \code{f} }
    \item{\code{name}:}{A string containing the name of the object}
    \item{\code{description}:}{A string containing a longer description
      of the object}
  }
}
\section{Methods}{
  \describe{
    \item{process(object, action, parameter)}{Apply the function
      represnted by \code{action} to the \code{Channel} object, updating
      the history appropriately. If the \code{parameter} is \code{NULL},
      then use the default value.
    }
    \item{summary(object, ...)}{Write out a summary of the object.}
  }
}
\value{
  The return value of the generic function \code{process} is always
  an object related to its \code{Channel} input, which keeps a record
  of its history. The precise class of the result depends on the
  function used to create the \code{Processor}.
}
\section{Pre-defined Processors}{
  The library comes with several \code{Processor} objects already
  defined; each one takes a \code{Channel} as input and produces a
  modified \code{Channel} as output.
  \describe{
    \item{\code{PROC.SUBTRACTOR}}{Subtracts a global constant (default:
      0) from the data vector in the \code{Channel}.}
    \item{\code{PROC.THRESHOLD}}{Truncates the data vector below,
      replacing the values below a threshold (default: 0) with the
      threshold value.}
    \item{\code{PROC.GLOBAL.NORMALIZATION}}{Normalizes the data vector
      in the \code{Channel} by dividing by a global constant.  If the
      parameter takes on its default value of 0, then divide by the 75th
      percentile.
    }
    \item{\code{PROC.LOG.TRANSFORM}}{Performs a log transformation of
      the data vector. The parameter specifies the base of the logarithm
      (default: 2).
    }
    \item{\code{PROC.MEDIAN.EXPRESSED.NORMALIZATION}}{Normalizes the
      data vector by dividing by the median of the expressed genes, where
      ``expressed'' is taken to mean ``greater than zero''.
    }
    \item{\code{PROC.SUBSET.NORMALIZATION}}{Normalizes the data vector
      by dividing by the median of a subset of genes.  When the
      parameter has a default value of 0, then this method uses the
      global median. Otherwise, the parameter should be set to a logical
      or numerical vector that selects the subset of genes to be used
      for normalization.
    }
    \item{\code{PROC.SUBSET.MEAN.NORMALIZATION}}{Normalizes the data vector
      by dividing by the mean of a subset of genes.  When the
      parameter has a default value of 0, then this method uses the
      global mean. Otherwise, the parameter should be set to a logical
      or numerical vector that selects the subset of genes to be used
      for normalization.
    }
  }
}
\author{Kevin R. Coombes <kcoombes@mdanderson.org>}
\seealso{ \code{\link{Channel}}, \code{\link{process}},
  \code{\link{Pipeline}}, \code{\link{CompleteChannel}}
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
nor <- process(subbed, PROC.GLOBAL.NORMALIZATION)
thr <- process(nor, PROC.THRESHOLD, 25)
processed <- process(thr, PROC.LOG.TRANSFORM, 2)

summary(processed)

par(mfrow=c(2,1))
plot(processed)
hist(processed)

par(mfrow=c(1,1))
image(processed)

rm(nor, thr, subbed, processed)
}
\keyword{manip}
