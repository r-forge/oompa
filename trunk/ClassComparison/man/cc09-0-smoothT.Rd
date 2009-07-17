\name{SmoothTtest-class}
\docType{class}
\alias{SmoothTtest}
\alias{SmoothTtest-class}
\alias{as.data.frame,SmoothTtest-method}
\alias{plot,SmoothTtest,missing-method}
\alias{summary,SmoothTtest-method}

\title{The SmoothTtest Class }
\description{
  Implements the smooth t-test for differential expression as developed
  by Baggerly and Coombes.
  }
\usage{
SmoothTtest(stats, aname = 'Group One', bname = 'Group Two',
 name = paste(aname, 'vs.', bname))
\S4method{as.data.frame}{SmoothTtest}(x, row.names=NULL, optional=FALSE)
\S4method{summary}{SmoothTtest}(object, ...)
\S4method{plot}{SmoothTtest,missing}(x, folddiff=3, goodflag=2, badch=4, ccl=0,
 name=x@name, pch='.', xlab='log intensity', ylab='log ratio', ...)
}
\arguments{
  \item{stats}{An object of the \code{\link{TwoGroupStats}} class.}
  \item{aname}{A character string; the name of the second group}
  \item{bname}{A character string; the name of the second group}
  \item{name}{A character string; the name of this object}
  \item{object}{A \code{SmoothTtest} object}
  \item{x}{A \code{SmoothTtest} object}
  \item{row.names}{See the base version of \code{\link{as.data.frame.default}}}
  \item{optional}{See the base version of
    \code{\link{as.data.frame.default}}}
  \item{folddiff}{A real number; the level of fold difference considered
    large enough to be indicated in the plots.}
  \item{goodflag}{A real number; the level (in standard deviation units)
    of the smooth t-statistic     considered large enough to be
    indicated in the plot.}
  \item{badch}{A real number; the level of variability in single groups
    considered large enough to be worrisome. See the \code{multiple}
    argument to the \code{plot} method in the \code{\link{SingleGroup}}
    class.} 
  \item{ccl}{A list containing objects of the
    \code{\link[PreProcess]{ColorCoding}} class. If left at its default
    value of zero, colors are chosen automatically.}
  \item{pch}{Default plotting character}
  \item{xlab}{Label for the x axis}
  \item{ylab}{Label for the y axis}
  \item{\dots}{The usual extra parameters to generic or plotting routines}
}
\details{
  In 2001 and 2002, Baggerly and Coombes developed the smooth t-test for
  finding differentially expressed genes in microarray data.  Along with
  many others, they began by log-transforming the data as a reasonable
  step in the direction of variance stabilization.  They observed,
  however, that the gene-by-gene standard deviations still seemed to
  vary in a systematic way as a function of the mean log intensity.  By
  borrowing strenght across genes and using \code{\link{loess}} to fit
  the observed standard deviations as a function of the mean, one
  presumably got a better estimate of the true standard deviation.

  These smooth estimates are computed for each of two groups of samples
  being compared.  They are then combined (gene-by-gene using the usual
  univariate formulas) to compute pooled "smooth" estimates of the standard
  deviation.  These smooth estimates are then used in gene-by-gene
  t-tests.

  The interesting question then arises of how to compute and interpret
  p-values associated to these individual tests.  The liberal
  argument asserts that, because smoothing uses data from hundreds
  of measurements to estiamte the standard deviation, it can effectively
  be treated as "known" in the t-tests, which should thus be compared
  against the normal distribution.  A conservative argument claims
  that the null distribution should still be the t-distribution with the
  degrees of freedom determined in the usual way by the number of
  samples.  The truth probably lies somewhere in between, and is
  probably best approximated by some kind of permutation test.  In this
  implementation, we take the coward's way out and don't provide any of
  those alternatives. You have to extract the t-statistics (from the
  \code{smooth.t.statistics} slot of the object) and compute your own
  p-values in your favorite way.  If you base the computations on a
  theoretical model rather than a permutation test, then the
  \code{\link{Bum}} class provides a convenient way to account for
  multiple testing.
}
\section{Creating Objects}{
  In practice, users will first use a data frame and a classification
  vector (or an \code{ExpressionSet}) to construct an object of the
  \code{\link{TwoGroupStats}} object.  This object can then be handed
  directly to the \code{SmoothTtest} function to perform the smooth
  t-test.
}
\section{Slots}{
  \describe{
    \item{\code{one}:}{An object of the \code{SingleGroup} class
      representing a loess smooth of standard deviation as a function of
      the mean in the first group of samples.}
    \item{\code{two}:}{An object of the \code{SingleGroup} class
      representing a loess smooth of standard deviation as a function of
      the mean in the second group of samples.}
    \item{\code{smooth.t.statistics}:}{A numeric vector containing the
      smooth t-statistics}
    \item{\code{fit}:}{A data.frame containing two columns, \code{x} and
      \code{y}, containing the smooth estimates of the pooled standard deviation}
    \item{\code{dif}:}{A numeric vector of the differences in mean
      values between the two groups}
    \item{\code{avg}:}{A numeric vector of the overall mean value}
    \item{\code{aname}:}{A character string; the name of the first group}
    \item{\code{bname}:}{A character string; the name of the second group}
    \item{\code{name}:}{A character string; the name of this object}
    \item{\code{stats}:}{The \code{TwoGroupStats} object that was used
      to create this object. }
  }
}
\section{Methods}{
  \describe{
    \item{as.data.frame(x, row.names=NULL, optional=FALSE)}{Convert the
      object into a data frame suitable for printing or exporting.}
    \item{summary(object, ...)}{Write out a summary of the object.}
    \item{plot(x, folddiff=3, goodflag=2, badch=4, ccl=0, name=x@name,
      pch='.', xlab='log intensity', ylab='log ratio', ...)}{Create a
      set of six plots.  The first two plots are the QC plots from the
      \code{\link{SingleGroup}} objects representing the two groups of
      samples.  The third plot is a scatter plot comparing the means in
      the two groups.  The fourth plot is Bland-Altman plot of the
      overall mean against the difference in means (also known
      colloquially as an M-vs-A plot).  The fifth plot is a histogram of
      the smooth t-statistics.  The final plot is a scatter plot of the
      smooth t-statistics as a function of the mean intensity. Colors in
      the plots are controlled by the curent values of
      \code{\link[PreProcess:semanticColors]{COLOR.BORING}},
      \code{\link[PreProcess:semanticColors]{COLOR.SIGNIFICANT}},
      \code{\link[PreProcess:semanticColors]{COLOR.BAD.REPLICATE}},
      \code{\link[PreProcess:semanticColors]{COLOR.WORST.REPLICATE}},
      \code{\link[PreProcess:semanticColors]{COLOR.FOLD.DIFFERENCE}},
      \code{\link[PreProcess:semanticColors]{COLOR.CENTRAL.LINE}}, AND
      \code{\link[PreProcess:semanticColors]{COLOR.CONFIDENCE.CURVE}}.
    }
  }
}
\references{
  Baggerly, K.A., Coombes, K.R., Hess, K.R., Stivers, D.N., Abruzzo,
  L.V., Zhang, W.  Identifying differentially expressed genes in cDNA
  microarray experiments.  J Comp Biol. 8:639-659, 2001.
  
  Coombes, K.R., Highsmith, W.E., Krogmann, T.A., Baggerly, K.A.,
  Stivers, D.N., Abruzzo, L.V.   Identifying and quantifying sources of
  variation in microarray data using high-density cDNA membrane arrays.
  J Comp Biol. 9:655-669, 2002.

  Altman DG, Bland JM. Measurement in Medicine: the Analysis of Method
  Comparison Studies. The Statistician, 1983; 32: 307-317. 
}
\author{Kevin R. Coombes <kcoombes@mdanderson.org>}
\seealso{
  \code{\link{Bum}}, \code{\link{MultiTtest}},
  \code{\link{SingleGroup}}, \code{\link{TwoGroupStats}}.
}
\examples{
bogus <- matrix(rnorm(30*1000, 8, 3), ncol=30, nrow=1000)
splitter <- rep(FALSE, 30)
splitter[16:30] <- TRUE
x <- TwoGroupStats(bogus, splitter)
y <- SmoothTtest(x)

opar <- par(mfrow=c(2, 3), pch='.')
plot(y, badch=2, goodflag=1)
par(opar)

# cleanup
rm(bogus, splitter, x, y, opar)
}
\keyword{classes}
\keyword{htest}
\keyword{multivariate}
