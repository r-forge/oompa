\name{Bum-class}
\docType{class}
\alias{Bum}
\alias{Bum-class}
\alias{countSignificant,Bum-method}
\alias{cutoffSignificant,Bum-method}
\alias{hist,Bum-method}
\alias{image,Bum-method}
\alias{selectSignificant,Bum-method}
\alias{summary,Bum-method}
\alias{likelihoodBum}

\title{The Bum Class}

\description{
  The \code{Bum} class is used to fit a beta-uniform mixture model to a
  set of p-values.
}
\usage{
Bum(pvals, ...)
\S4method{summary}{Bum}(object, tau=0.01, ...)
\S4method{hist}{Bum}(x, res=100, xlab='P Values', main='', ...)
\S4method{image}{Bum}(x, ...)
\S4method{cutoffSignificant}{Bum}(object, alpha, by='FDR', ...)
\S4method{selectSignificant}{Bum}(object, alpha, by='FDR', ...)
\S4method{countSignificant}{Bum}(object, alpha, by='FDR', ...)
likelihoodBum(object)
}
\arguments{
  \item{pvals}{A numeric vector containing values between 0 and 1}
  \item{object}{A \code{Bum} object}
  \item{tau}{A real number between 0 and 1, representing a cutoff on the
    p-values.}
  \item{x}{A \code{Bum} object}
  \item{res}{A positive integer; the resolution at which to plot the
    fitted distribution curve.}
  \item{xlab}{Label for the x axis}
  \item{main}{Graph title}
  \item{alpha}{Either the false discovery rate (if \code{by = 'FDR'}) or
  the posterior probability (if \code{by = 'EmpiricalBayes'})}
  \item{by}{String denoting the method to use for determining
    cutoffs. The chioces are 'FDR', 'FalseDiscovery', or
    'EmpiricalBayes'. Since the test is implemented with
    \code{\link{match.arg}}, unique abbreviations also work.}
  \item{\dots}{All methods are defined to accept additional arguments in
    order to allow flexibility in designing derived classes.  The usual
    graphical parameters can be supplied to \code{hist} and \code{image}.}
}
\section{Creating Objects}{
  Although objects can be created directly using \code{new}, the most
  common usage will be to pass a vector of p-values to the
  \code{Bum} function. 
}
\section{Slots}{
  \describe{
    \item{\code{pvals}:}{The vector of p-values used to construct the object.}
    \item{\code{ahat}:}{Model parameter}
    \item{\code{lhat}:}{Model parameter}
    \item{\code{pihat}:}{Model parameter}
  }
}
\section{Methods}{
  \describe{
    \item{summary(object, tau=0.01, ...)}{For each value of the p-value
      cutoff \code{tau}, computes estimates of the fraction of true
      positives (TP), false negatives (FN), false positives (FP), and
      true negatives (TN).}
    \item{hist(x, res=100, xlab='P Values', main='', ...)}{ Plots a
      histogram of the object, and overlays (1) a straight line to indicate
      the contribution of the uniform component and (2) the fitted
      beta-uniform distribution from the observed values. Colors in the
      plot are controlled by
      \code{\link[PreProcess:semanticColors]{COLOR.EXPECTED}} and 
      \code{\link[PreProcess:semanticColors]{COLOR.OBSERVED}}.
    }
    \item{image(x, ...)}{Produces four plots in a 2x2 layout: (1) the
      histogram produced by \code{hist}; (2) a plot of cutoffs against
      the desired false discovery rate; (3) a plot of cutoffs against
      the posterior probability of coming from the beta component; and
      (4) an ROC curve.}
    \item{cutoffSignificant(object, alpha, by='FDR', ...)}{Computes the
      cutoff needed for significance, which in this case means arising
      from the beta component rarther than the uniform component of the
      mixture. Significance is specified either by the false discovery
      rate (when \code{by = 'FDR'} or \code{by = 'FalseDiscovery'}) or
      by the posterior probability (when \code{by = 'EmpiricalBayes'})}
    \item{selectSignificant(object, alpha, by='FDR', ...)}{Uses
      \code{cutoffSignificant} to determine a logical vector that
      indicates which of the p-values are significant.}
    \item{countSignificant(object, alpha, by='FDR', ...)}{Uses
      \code{selectSignificant} to count the number of significant
      p-values.}
  }
}
\value{
  Graphical functions (\code{hist} and \code{image}) invisibly return the
  object on which they were invoked.

  The \code{cutoffSignficant} method returns a real number between zero
  and one. P-values below this cutoff are considered statistically
  significant at either the specified false discovery rate or at the
  specified posterior probability.

  The \code{selectSignficant} method returns a vector of logical values
  whose length is equal to the length of the vector of p-values that was
  used to construct the \code{Bum} object.  True values in the return
  vector mark the statistically signficant p-values.

  The \code{countSignificant} method returns an integer, the number of
  statistically significant p-values.

  The \code{summary} method returns an object of class
  \code{\link{BumSummary}}.
}
\details{
  The BUM method was introduced by Stan Pounds and Steve Morris,
  although it was simultaneously discovered by several other
  researchers.  It is generally applicable to any analysis of microarray
  or proteomics data that performs a separate statistical hypothesis
  test for each gene or protein, where each test produces a p-value that
  would be valid if the analyst were only performing one statistical
  test.  When performing thousands of statistical tests, however, those
  p-values no longer have the same interpretation as Type I error
  rates.  The idea behind BUM is that, under the null hypothesis that
  none of the genes or proteins is interesting, the expected
  distribution of the set of p-values is uniform.  By contrast, if some
  of the genes are interesting, then we should see an overabundance of
  small p-values (or a spike in the histogram near zero).  We can model
  the alternative hypothesis with a beta distribution, and view the set
  of all p-values as a mixture distribution.

  Fitting the BUM model is straightforward, using a nonlinear optimizer
  to compute the maximum likelihood parameters.  After the model has
  been fit, one can easily determine cutoffs on the p-values that
  correspond to desired false discovery rates.  Alternatively, the
  original Pounds and Morris paper shows that their results can be
  reinterpreted to recover the empirical Bayes method introduced by
  Efron and Tibshirani.  Thus, one can also determine cutoffs by
  specifying a desired posterior probability of signficance.
}
\references{
  Pounds S, Morris SW. Estimating the occurrence of false
  positives and false negatives in microarray studies by approximating
  and partitioning the empirical distribution of
  p-values. Bioinformatics. 2003 Jul 1;19(10):1236-42. 

  Benjamini Y, Hochberg Y.  Controlling the false discovery rate: a
  practical and powerful approach to multiple testing. J Roy Statist Soc
  B, 1995; 57: 289-300.

  Efron B, Tibshirani R: Empirical bayes methods and false discovery
  rates for microarrays. Genet Epidemiol 2002, 23: 70-86.
}
\author{Kevin R. Coombes <kcoombes@mdanderson.org>}
\seealso{Two classes that produce lists of p-values that can (and often
  should) be analyzed using BUM are \code{\link{MultiTtest}} and
  \code{\link{MultiLinearModel}}. Also see \code{\link{BumSummary}}. 
}
\examples{
fake.data <- c(runif(700), rbeta(300, 0.3, 1))
a <- Bum(fake.data)
hist(a, res=200)

alpha <- (1:25)/100
plot(alpha, cutoffSignificant(a, alpha, by='FDR'),
     xlab='Desired False Discovery Rate', type='l',
     main='FDR Control', ylab='Significant P Value')

GAMMA <- 5*(10:19)/100
plot(GAMMA, cutoffSignificant(a, GAMMA, by='EmpiricalBayes'),
     ylab='Significant P Value', type='l',
     main='Empirical Bayes', xlab='Posterior Probability')

b <- summary(a, (0:100)/100)
be <- b@estimates
sens <- be$TP/(be$TP+be$FN)
spec <- be$TN/(be$TN+be$FP)
plot(1-spec, sens, type='l', xlim=c(0,1), ylim=c(0,1), main='ROC Curve')
points(1-spec, sens)
abline(0,1)

image(a)

countSignificant(a, 0.05, by='FDR')
countSignificant(a, 0.99, by='Emp')

#cleanup
rm(a, b, be, alpha, GAMMA, sens, spec, fake.data)
}
\keyword{multivariate}
\keyword{htest}
\keyword{classes}

