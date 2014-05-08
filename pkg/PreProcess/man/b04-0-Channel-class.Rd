\name{Channel-class}
\alias{Channel-class}
\alias{Channel}
\alias{plot,Channel,missing-method}
\alias{hist,Channel-method}
\alias{summary,Channel-method}
\alias{print,Channel-method}
\alias{show,Channel-method}
\alias{image,Channel-method}
\docType{class}
\title{Class "Channel"}
\description{
  An object of the \code{Channel} class represents a single kind of
  measurement performed at all spots of a microarray channel. These
  objects are essentially just vectors of data, with length equal to the
  number of spots on the microarray, with some extra metadata attached.
}
\usage{
Channel(parent, name, type, vec)
\S4method{plot}{Channel,missing}(x, y, \dots)
\S4method{hist}{Channel}(x, breaks=67, xlab=x@name, main=x@parent, \dots)
\S4method{summary}{Channel}(object, \dots)
\S4method{print}{Channel}(x, \dots)
\S4method{show}{Channel}(object)
\S4method{image}{Channel}(x, main=x@name, sub=NULL, \dots)
}
\arguments{
  \item{parent}{character string representing the name of a parent object from
    which this object was derived}
  \item{name}{character string with a displayable name for this object}
  \item{type}{object of class \code{\link{ChannelType}}}
  \item{vec}{numeric vector}
  \item{x}{object of class \code{Channel}}
  \item{y}{nothing; the new Rd format requires documenting missing parameters}
  \item{breaks}{see the documentation for the default \code{\link[graphics]{hist}}}
  \item{xlab}{character string specifying the label for x axis}
  \item{main}{character string specifying the main title for the plot}
  \item{sub}{character string specifying subtitle for the plot}
  \item{object}{object of class \code{Channel}}
  \item{\dots}{extra arguments for generic or plotting routines}
}
\details{
  As described in the help pages for \code{\link{ChannelType}}, each
  microarray hybridization experiment produces one or more channels of
  data. \code{Channel} objects represent a single measurement performed
  at spots in one microarray channel. The raw data from a full experiment
  typically contains multiple measurements in multiple channels.

  The full set of measurements is often highly processed (by, for example,
  background subtraction, normalization, log transformation, etc.) before it
  becomes useful. We have added a \code{history} slot that keeps track of how
  a \code{Channel} was produced. By allowing each object to maintain a record
  of its history, it becomes easier to document the processing when writing up
  the methods for reports or papers. The \code{history} slot of the object is
  updated using the generic function \code{\link{process}} together with a
  \code{\link{Processor}} object.
}
\section{Slots}{
  \describe{
    \item{\code{parent}:}{character string representing the name of a parent
      object from which this object was derived.
    }
    \item{\code{name}:}{character string with a displayable name for this
      object}
    \item{\code{type}:}{object of class \code{\link{ChannelType}}}
    \item{\code{x}:}{numeric vector}
    \item{\code{history}:}{list that keeps a record of the calls used
      to produce this object}
  }
}
\section{Methods}{
  \describe{
    \item{print(object, \dots)}{Print all the data on the object. Since
      this includes the entire data vector, you rarely want to do this.}
    \item{show(object)}{Print all the data on the object. Since
      this includes the entire data vector, you rarely want to do this.}
    \item{summary(object, \dots)}{Write out a summary of the object.}
    \item{plot(object, \dots)}{Produce a scatter plot of the measurement
      values in the slot \code{x} of the \code{object} against their
      index , which serves as a surrogate for the position on the
      microarray. Additional graphical parameters are passed along.}
    \item{hist(object, \dots)}{Produce a histogram of the data values in slot
      \code{x} of the \code{object}. Additional graphical parameters are
      passed along.}
    \item{image(object, \dots)}{This method produces a two-dimensional
      "cartoon" image of the measurement values, with the position in
      the cartoon corresponding to the two-dimensional arrangement of
      spots on the actual microarray. Additional graphical parameters
      are passed along.}
  }
}
\value{
  The \code{print}, \code{hist}, and \code{image} methods all invisibly
  return the \code{Channel} object on which they were invoked.

  The \code{print} and \code{summary} methods return nothing.
}
\author{
  Kevin R. Coombes \email{krc@silicovore.com},
  P. Roebuck \email{proebuck@mdanderson.org}
}
\seealso{
  \code{\link{ChannelType}},
  \code{\link{process}},
  \code{\link{Processor}}
}
\examples{
showClass("Channel")

## simulate a moderately realistic looking microarray
nc <- 100			# number of rows
nr <- 100			# number of columns
v <- rexp(nc*nr, 1/1000)	# "true" signal intensity (vol)
b <- rnorm(nc*nr, 80, 10)	# background noise
s <- sapply(v-b, max, 1)	# corrected signal intensity (svol)
ct <- ChannelType('user', 'random', nc, nr,  'fake')
raw <- Channel(name='fraud', type=ct, parent='', vec=v)
subbed <- Channel(name='fraud', parent='', type=ct, vec=s)
rm(nc, nr, v, b, s)		# clean some stuff

summary(subbed)
summary(raw)

par(mfrow=c(2,1))
plot(raw)
hist(raw)

par(mfrow=c(1,1))
image(raw)

## finish the cleanup
rm(ct, raw, subbed)
}
\keyword{classes}
\keyword{manip}

