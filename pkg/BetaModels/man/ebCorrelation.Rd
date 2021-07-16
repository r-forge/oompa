\name{ebCorrelation-class}
\alias{ebCorrelation-class}
\alias{ebCorrelation}
\alias{hist,ebCorrelation-method}
\alias{plot,ebCorrelation,missing-method}
\alias{cutoffSignificant,ebCorrelation-method}
\alias{selectSignificant,ebCorrelation-method}
\alias{countSignificant,ebCorrelation-method}
\alias{summary,ebCorrelation-method}


\title{Class "ebCorrelation"}
\description{Fit an empirical Bayes model to detemine which values in a
  large collection of correlation coefficients are significant.
}
\usage{
ebCorrelation(ss, nObs, nPoints = 500)
\S4method{hist}{ebCorrelation}(x,
           xlab='Correlation', ylab='Prob(Different | Y)', main='',
           highlight='purple', lowlight='blue', \dots)
\S4method{plot}{ebCorrelation,missing}(x,
                  prior=1, significance=0.9, ylim=c(-0.5, 1),
                  xlab='Correlation', ylab='Prob(Unusual | Rho)',
                  highlight='purple', \dots)
\S4method{summary}{ebCorrelation}(object, prior=1, significance=0.9, \dots)
\S4method{cutoffSignificant}{ebCorrelation}(object, prior, significance, \dots)
\S4method{selectSignificant}{ebCorrelation}(object, prior, significance, \dots)
\S4method{countSignificant}{ebCorrelation}(object, prior, significance, \dots)
}
\arguments{
  \item{ss}{A numerical vector containing correlation coefficinets
    between -1 and 1.}
  \item{nObs}{A numerical vector of length one, the number of objects used in every
    computation of correlation coefficients.}
  \item{nPoints}{the number of points at which to estimate the distribution.}
  \item{object}{object of class \code{ebCorrelation}}
  \item{x}{object of class \code{ebCorrelation}}
  \item{xlab}{Graphical parameter.}
  \item{ylab}{Graphical parameter.}
  \item{ylim}{Graphical parameter.}
  \item{main}{Graphical parameter.}
  \item{\dots}{Optional extra parameters, either graphical or for the
    significance functions.}
  \item{prior}{A real number between 0 and 1; the prior probability that
    a correlation coefficinet is not significant.}
  \item{significance}{A real number between 0 and 1; the posterior
    probability betyond which a correlation coefficinet will be called
    significant.}
  \item{highlight}{Character string denoting a color.}
  \item{lowlight}{Character string denoting a color.}
}
\section{Creating Objects}{
  Although objects can be created directly using \code{new}, the most
  common usage will be to pass a vector of correlation coefficients to the
  \code{ebCorrelation} function.
}
\section{Slots}{
  \describe{
    \item{\code{correlation}:}{numeric vector of correlation coefficients.}
    \item{\code{nObservation}:}{the number of sample observations used
      to compute correlations.}
    \item{\code{xvals}:}{vector of the x-axis grid points at which to compute the
      posterior probability; see Details.}
    \item{\code{pdf}:}{vector of the ermpirically estimated probability densities
      at \code{xvals}; see Details.}
    \item{\code{theoretical.pdf}:}{vector of the theoretical probability densities
      at \code{xvals}; see Details.}
    \item{\code{unravel}:}{Matrix of posterior probabilities.}
    \item{\code{call}:}{A \code{call} object recording how the
      constructior function was invoked.} 
  }
}
\section{Methods}{
  \describe{
    \item{summary(object, \dots)}{Prints a summary of the ebCorrelation
      object.  This includes (1) the maximum a posterior coordinates on
      x-y-space, (2) the usual alpha-beta parameters for the Beta
      distribution, and (3) the mean and variance.}
    \item{image(x, col, \dots)}{Plots an ikmage of the posterior
      probabilities using te specified color map.  The point with the
      maximum posterior probability is marked in red.}
  }
}
\value{
  The \code{ebCorrelation} constructor returns an object of the indicated
  class.

  TBD
}
\details{
  TBD
}
\references{
  Efron's paper on empirical Bayes for differential exprewssion.
}
\author{
  Kevin R. Coombes \email{krc@silicovore.com}
}
\examples{
showClass("ebCorrelation")
set.seed(12345)
cc <- c(rbeta(4600, 24, 24), rbeta(400, 8, 8))
rr <- 2*cc-1
fit <- ebCorrelation(rr, 51)
hist(fit)
plot(fit, prior = 0.85)
countSignificant(fit, prior = 0.85, significance = 0.8)
cutoffSignificant(fit, prior = 0.85, significance = 0.8)
summary(fit)
summary(fit, prior = 0.85, significance = 0.8)

}
\keyword{classes}
\keyword{htest}
