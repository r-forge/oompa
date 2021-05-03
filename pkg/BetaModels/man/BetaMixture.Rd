\name{BetaMixture-class}
\alias{BetaMixture-class}
\alias{BetaMixture}
\alias{summary,BetaMixture-method}
\alias{hist,BetaMixture-method}

\title{Class "BetaMixture"}
\description{Model data on the interval [0,1] as a mixture of an
  arbitrary number of beta distributions.
}
\usage{
BetaMixture(datavec, K = 2, forever = 100, epsilon = 0.001, debug = FALSE)
\S4method{summary}{BetaMixture}(object, \dots)
\S4method{hist}{BetaMixture}(x, mixcols = 1:7, \dots)
}
\arguments{
  \item{datavec}{The observed vector of data points (between 0 and 1).}
  \item{K}{The number of mixture components to fit.}
  \item{forever}{The maximum number of iterations of the algorithm.}
  \item{epsilon}{The minimum change of parameters that must be observed
    in order to continue iterating.}
  \item{debug}{A local valuecontrolling whethe to print out intermediate steps.}
  \item{x}{In the \code{hist} method, an object of class \code{BetaMixture}.}
  \item{mixcols}{A vector of colors used to show differnt
    components. Recycled as necessary.}
  \item{object}{Object of class \code{BetaMixture}.}
  \item{\dots}{Extra arguments for generic routines.}
}
\section{Creating Objects}{
  Although objects can be created directly using \code{new}, the most
  common usage will be to pass a vector of unit interval values to the
  \code{BetaMixture} function.
}
\section{Slots}{
  \describe{
    \item{\code{datavec}:}{input data vector of values between 0 and 1.}
    \item{\code{phi}:}{vector of non-negative component frequencies;  will sum to 1.}
    \item{\code{mle}:}{\eqn{$2\times K$}{2xK} matrix of parameters of the $K$ beta
      mixture components.}
    \item{\code{Z}:}{\eqn{$N\times K$}{NxK} matrix of the probability for each of
      the $N$ observations to arise from beta mixture component $K$. Each
      row will sum to 1.}
    \item{\code{negloglike}:}{The negative log likelihood of the fitted model.}
    \item{\code{converged}:}{A logical value indicating where the EM
      algorithm converged.}
  }
}
\section{Methods}{
  \describe{
    \item{summary(object, \dots)}{Prints a summary of the \code{BetaMixture}
      object.  This includes (1) the frequency/weight of each component, (2)
      the parameters of each beta distribution, (3) an indicator of
      convergence, and (4) the log-likelihood of the fimal model.}
  }
}
\value{
  The \code{BetaMixture} constructor returns an object of the indicated
  class.

  The graphical method \code{hist} and the \code{summary} method
  invisibly return the object on which it was invoked.
}
\details{
  Given a data set consisting of values in the interval [0,1], we want to
  model it as a mixture of \eqn{$K$}{K} beta distributions:
  \deqn{f(x) = \sum_{i=1}^K \phi_i B(\nu_i, \omega_i)}
  {f(X) = sum phi[u] B(nu[i], omega[i])}
  where the non-negative frequencies or weights of the components sum to
  1 (that is, \eqn{\sum_{i=1}^K \phi_i = 1}{sum phi[i] = 1}).

  We fit such a model using en Expectation-Maximation (EM) algorithm.
  To accomplish this goal, we must introduce a matrix \eqn{$Z$}{Z} of
  latent assignments of each of the $N$ observation to the $K$
  components. In the current implementation, we intialize teh
  assignments by dividing the range of the observed data into $K$
  subintervals of equal width, and assign the elemnts in the Kth
  subinterval to component $K$. (\strong{Warning}: This may change in later
  versions.) We then alternate between the M-step (finding a maximum
  likelihood estimate of the parameters of the beta distributions given
  Z), and the E-step (estimating the weights \eqn{\phi}{phi} by
  averaging the columns of Z). Iterations continue up to the number
  specified by \code{forever} times, or until the log likelihood changes
  by less than the amount specified by \code{epsilon}.
}
\author{
  Kevin R. Coombes \email{krc@silicovore.com}
}
\examples{
showClass("BetaMixture")
set.seed(73892)
datavec <- c(rbeta(130, 1, 4),
             rbeta(170, 7, 4))
model <- BetaMixture(datavec)
summary(model)
hist(model, breaks=35)
}
\keyword{classes}
\keyword{htest}
