\name{BetaRates-class}
\alias{BetaRates-class}
\alias{BetaRates}
\alias{summary,BetaRates-method}
\alias{image,BetaRates-method}
\alias{guessCenter}
\alias{samplePosteriorRates}

\title{Class "BetaRates"}
\description{Test for different proportions (rates) in different groups
  using a Bayesian model in which all rate parameters follow a beta
  distribution and are selected from a common hyperdistribution. 
}
\usage{
BetaRates(k, n, x=seq(-3,3, length=100), y=x)
\S4method{summary}{BetaRates}(object, \dots)
\S4method{image}{BetaRates}(x, col=greyscale(128), \dots)
samplePosteriorRates(br, nsamp=2000)
guessCenter(v)
}
\arguments{
  \item{object}{object of class \code{BetaRates}}
  \item{br}{object of class \code{BetaRates}}
  \item{x}{In the \code{image} method, an object of class
    \code{BetaRates}. In the \code{BetaRates} constructor, a vector of
    the x-axis grid points at which to compute the posterior
    probability; see Details.}
  \item{y}{vector of the y-axis grid points at which to compute the
    posterior probability; see Details.}
  \item{k}{vector of "success" counts}
  \item{n}{vector of all counts}
  \item{col}{vector containing the color map to use for the image}
  \item{nsamp}{Number of posterior samples to take}
  \item{v}{Vector of observed rates}
  \item{\dots}{extra arguments for generic or plotting routines}
}
\section{Creating Objects}{
  Although objects can be created directly using \code{new}, the most
  common usage will be to pass a vector of p-values to the
  \code{BetaRates} function. 
}
\section{Slots}{
  \describe{
    \item{\code{k}:}{vector of "success" counts.}
    \item{\code{n}:}{vector of all counts.}
    \item{\code{x}:}{vector of the x-axis grid points at which to compute the
    posterior probability; see Details.}
    \item{\code{y}:}{vector of the y-axis grid points at which to compute the
    posterior probability; see Details.}
    \item{\code{results}:}{Matrix of posterior probabilities.}
    \item{\code{logresults}:}{Matrix of log-transformed posterior
      probabilities.}
  }
}
\section{Methods}{
  \describe{
    \item{summary(object, \dots)}{Prints a summary of the betaRates
      object.  This includes (1) the maximum a posterior coordinates on
      x-y-space, (2) the usual alpha-beta parameters for the Beta
      distribution, and (3) the mean and variance.}
    \item{image(x, col, \dots)}{Plots an ikmage of the posterior
      probabilities using te specified color map.  The point with the
      maximum posterior probability is marked in red.}
  }
}
\value{
  The \code{BetaRates} constructor returns an object of the indicated
  class.

  The graphical method \code{image}) invisibly returns the object on
  which it was invoked.

  The \code{summary} method returns a vector with the maximum a
  posteriori parameters of the beta distribution.

  The \code{samplePosteriorRates} function returns a list with two
  components. The first component, \code{xy}, is an \code{nsamp}-by-2
  matrix with x-y values samples from the posterior distribution. The
  second component, \code{theta}, is an \code{nsamp}-by-\code{length(k)}
  matrix with posterior samples of the rates associated with each
  experiment supplied to the constructor.
  
  The \code{guessCenter} function returns a list with both x-y and
  alpha-beta coordinates of the naive (frequentist) estimate fo mthe
  overall Beta distribution parameters.
}
\details{
  TBD
}
\references{
  Gelman A, Carlin JB, Stern HS, Rubin DB. Bayesian Data Analysis,
  second edition. Chapman and Hall/CRC, Boca Raton, 2004. Section 5.3,
  pages 15-131.
}
\author{
  Kevin R. Coombes \email{krc@silicovore.com}
}
\examples{
showClass("BetaRates")
event <- c( 37,  4,  6,  1,  2, 10,  1, 13,   7,  1,  10)
total <- c(137, 18, 18, 26, 24, 45, 12, 43, 162, 78, 280)
guessCenter(event/total)
br <- BetaRates(event, total, x=seq(-3, 0, length=100), y=seq(0, 3, length=100))
image(br)
summary(br)
}
\keyword{classes}
\keyword{htest}
