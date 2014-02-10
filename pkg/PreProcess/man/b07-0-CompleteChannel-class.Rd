\name{CompleteChannel-class}
\alias{CompleteChannel-class}
\alias{CompleteChannel}
\alias{print,CompleteChannel-method}
\alias{summary,CompleteChannel-method}
\alias{channelize,CompleteChannel-method}
\alias{plot,CompleteChannel,missing-method}
\alias{image,CompleteChannel-method}
\alias{as.data.frame,CompleteChannel-method}
\alias{analyze,CompleteChannel-method}
\alias{process,CompleteChannel,Processor-method}
\alias{PROC.BACKGROUND}
\alias{PROC.SIGNAL}
\alias{PROC.CORRECTED.SIGNAL}
\alias{PROC.NEG.CORRECTED.SIGNAL}
\alias{PROC.SD.SIGNAL}
\alias{PROC.SIGNAL.TO.NOISE}
\docType{class}
\title{Class "CompleteChannel"}
\description{
  An object of the \code{CompleteChannel} class represents one channel (red or
  green) of a two-color fluorescence microarray experiment. Alternatively,
  it can also represent the entirety of a radioactive microarray experiment.
  Affymetrix experiments produce data with a somewhat different structure
  because they use multiple probes for each target gene.
}
\usage{
CompleteChannel(name, type, data)
\S4method{print}{CompleteChannel}(x, \dots)
\S4method{summary}{CompleteChannel}(object, \dots)
\S4method{as.data.frame}{CompleteChannel}(x, row.names=NULL, optional=FALSE)
\S4method{plot}{CompleteChannel,missing}(x, main=x@name, useLog=FALSE, \dots)
\S4method{image}{CompleteChannel}(x, \dots)
\S4method{analyze}{CompleteChannel}(object, useLog=FALSE, \dots)
\S4method{process}{CompleteChannel,Processor}(object, action, parameter)
\S4method{channelize}{CompleteChannel}(object, \dots)
}
\arguments{
  \item{name}{character string specifying the name of the object}
  \item{type}{object of class \code{\link{ChannelType}}}
  \item{data}{data frame. For the pre-defined ``extraction''
    processors to work correctly, this should include columns called
    \code{vol}, \code{bkgd}, \code{svol}, \code{SD}, and \code{SN}.}
  \item{x}{object of class \code{CompleteChannel}}
  \item{object}{object of class \code{CompleteChannel}}
  \item{main}{character string specifying the title for the plot}
  \item{useLog}{logical scalar. If \code{TRUE}, convert to logarithmic values.}
  \item{action}{object of class \code{Processor} used to process a
    \code{CompleteChannel}}
  \item{parameter}{any object that makes sense as a parameter to the
    function represented by the \code{Processor} \code{action}}
  \item{row.names}{See \code{\link{as.data.frame}}}
  \item{optional}{See \code{\link{as.data.frame}}}
  \item{\dots}{extra arguments for generic or plotting routines}
}
\details{
  The names come from the default column names in the ArrayVision software
  package used at M.D. Anderson for quantifying glass or nylon microarrays.
  Column names used by other software packages should be mapped to these.
}
\section{Slots}{
  \describe{
    \item{\code{name}:}{character string containing the name of the object}
    \item{\code{type}:}{object of class \code{\link{ChannelType}}}
    \item{\code{data}:}{data frame}
    \item{\code{history}:}{list that keeps a record of the calls used
      to produce this object}
  }
}
\section{Methods}{
  \describe{
    \item{print(x, \dots)}{Print all the data on the object. Since
      this includes the data frame, you rarely want to do this.}
    \item{summary(object, \dots)}{Write out a summary of the object.}
    \item{as.data.frame(x,row.names=NULL, optional=FALSE)}{Convert the
      \code{CompleteChannel} object into a data frame. As you might
      expect, this simply returns the data frame in the \code{data} slot
      of the object.}
    \item{plot(x, useLog=FALSE, \dots)}{Produces three estimated density
      plots: one for the signal, one for the background, and one for
      the background-corrected signal. Additional graphical parameters
      are passed along. The logical flag \code{useLog} determines
      whether the data are log-transformed before estimating and
      plotting densities.} 
    \item{analyze(object, useLog=FALSE, \dots)}{This method computes the
      estimated probability density functions for the three data
      components (signal, background, and background-corrected signal),
      and returns them as a list. }
    \item{image(object, \dots)}{Uses the image method for
      \code{\link{Channel}} objects to produce geographically aligned
      images of the log-transformed intensity and background estimates.}
    \item{channelize(object, \dots)}{character string giving the name of the
      class of a channel that is produced when you process a
      \code{CompleteChannel} object.}
    \item{process(object, action, parameter=NULL)}{Use the
      \code{Processor} \code{action} to process the
      \code{CompleteChannel} \code{object}. Returns an object of the
      class described by \code{channelize}, which defaults to
      \code{Channel}.}
  }
}
\value{
  The \code{analyze} method returns a list of three density functions.

  The return value of the \code{process} function depends on the
  \code{\link{Processor}} performing the action, but is typically a
  \code{\link{Channel}} object.

  Graphical methods invisibly return the object on which they were
  invoked. 
}
\section{Pre-defined Processors}{
  The library comes with several \code{Processor} objects already
  defined; each one takes a \code{CompleteChannel} as input, extracts a
  single value per spot, and produces a \code{Channel} as output.
  \describe{
    \item{\code{PROC.BACKGROUND}}{Extract the vector of local
      background measurements.}
    \item{\code{PROC.SIGNAL}}{Extract the vector of foreground signal
      intensity measurements. }
    \item{\code{PROC.CORRECTED.SIGNAL}}{Extract the vector of
      background-corrected signal measurements. Note that many software
      packages automatically truncate these value below at zero, so this
      need not be the same as SIGNAL - BACKGROUND.}
    \item{\code{PROC.NEG.CORRECTED.SIGNAL}}{Extract the vector of
      background-corrected signal intensities by subtracting the local
      background from the observed foreground, without truncation.}
    \item{\code{PROC.SD.SIGNAL}}{Extract the vector of pixel standard
      deviations of the signal intensity.}
    \item{\code{PROC.SIGNAL.TO.NOISE}}{Extract the vector of
      signal-to-noise ratios, defined as CORRECTED.SIGNAL divided by the
      standard deviation of the background pixels.}
  }
}
\author{
  Kevin R. Coombes \email{krc@silicovore.com},
  P. Roebuck \email{proebuck@mdanderson.org}
}
\seealso{
  \code{\link{process}},
  \code{\link{Processor}},
  \code{\link{Pipeline}},
  \code{\link{Channel}},
  \code{\link{as.data.frame}}
}
\examples{
showClass("CompleteChannel")

## simulate a complete channel object
v <- rexp(10000, 1/1000)
b <- rnorm(10000, 60, 6)
s <- sapply(v-b, function(x) {max(0, x)})
ct <- ChannelType('user', 'random', 100, 100, 'fake')
x <- CompleteChannel(name='fraud', type=ct,
                      data=data.frame(vol=v, bkgd=b, svol=s))
rm(v, b, s, ct)

summary(x)

opar <- par(mfrow=c(2,3))
plot(x)
plot(x, main='Log Scale', useLog=TRUE)
par(opar)

opar <- par(mfrow=c(2,1))
image(x)
par(opar)

b <- process(x, PROC.NEG.CORRECTED.SIGNAL)
summary(b)

q <- process(b, PIPELINE.STANDARD)
summary(q)

q <- process(x, PIPELINE.MDACC.DEFAULT)
summary(q)

## cleanup
rm(x, b, q, opar)
}
\keyword{classes}
\keyword{manip}

