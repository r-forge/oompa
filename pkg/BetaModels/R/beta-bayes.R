# use empirical bayes to determine significance of correlation coefficients

setClass("ebCorrelation",
         representation(correlations='numeric',
                        nObservations='numeric',
                        xvals='numeric',
                        pdf='numeric',
                        theoretical.pdf='numeric',
                        unravel='numeric',
                        call='call'))

.probUnusual <- function(object, p0) {
  1-p0*object@theoretical.pdf/object@unravel
}

# input:
#  ss = is a vector of computed correlation coefficients
# nObs = number of observations used for each cor coef
#  nPoints = number of points at which to fit the distribution
ebCorrelation <- function(ss, nObs, nPoints=500) {
  call <- match.call()
  # start by estimating the empirical distribution
  eps <- 1/(2*nPoints)
  xvals <- seq(-1-eps, 1+eps, length=nPoints)
  pdf <- hist(ss, breaks=xvals, prob=TRUE, plot=FALSE)$density
  # get the theoretical distribution from a beta(M,M)
  M <- (nObs-3)/2
  # have to evaluate at the midpoint of the intervals
  xvals <- xvals[2:length(xvals)]-diff(xvals[1:2])
  # need to shift the interval back and correct using Jacobian
  theo <- dbeta((xvals+1)/2, M, M)/2
  # now we look at differnce between empirical and theoretical after log transofrm
  Y <- log(pdf)-log(theo)
  click <- pdf != 0 & theo != 0
  Y <- Y[click]
  X <- xvals[click]
  # fit this with a spline
  Z <- lm(Y ~ bs(X, df=25))
  YP <- predict(Z, data.frame(X=xvals))
  # convert back to the original scale
  unravel <- exp(YP+log(theo))
  new("ebCorrelation", correlations=ss, nObservations=nObs,
      xvals=xvals, pdf=pdf, theoretical.pdf=theo, unravel=unravel,
      call=call)
}


setMethod('hist', 'ebCorrelation', function(x,
           xlab='Correlation',
           ylab='Prob(Different | Y)', main='',
           highlight='purple', lowlight='blue', ...) {
  top <- max(c(x@unravel, x@theoretical.pdf))
  hist(x@correlations, probability=TRUE, breaks=100, ylim=c(0, top),
       xlim=c(min(x@xvals), max(x@xvals)), xlab=xlab, main=main)
  lines(x@xvals, x@theoretical.pdf, col=lowlight, lwd=2)
  lines(x@xvals, x@unravel, col=highlight, lwd=2)
  legend(min(x@xvals), max(x@unravel), c('Empirical', 'Theoretical'),
         col=c(highlight, lowlight), lwd=2)
  invisible(x)
})

setMethod('plot', signature('ebCorrelation', 'missing'), function(x,
                  prior=1, significance=0.9, ylim=c(-0.5, 1),
                  xlab='Correlation',
                  ylab='Prob(Unusual | Rho)',
                  highlight='purple', ...) {
  plot(c(min(x@xvals), max(x@xvals)), c(-0.5, 1), type='n',
       xlab=xlab, ylab=ylab, ylim=ylim, ...)
  toss <- unlist(lapply(prior, function(p, o) {
    lines(o@xvals, .probUnusual(o, p), err=-1)
  }, x))
  abline(h=significance, col=highlight)
  abline(h=0)
  invisible(x)
})

setMethod('cutoffSignificant', signature('ebCorrelation'),
          function(object, prior, significance, ...) {
            z <- .probUnusual(object, prior)
            x <- object@xvals[z < significance]
            list(low=min(x, na.rm=TRUE), high=max(x, na.rm=TRUE))
          })

setMethod('selectSignificant', signature('ebCorrelation'),
          function(object, prior, significance, ...) {
            stats <- object@correlations
            lh <- cutoffSignificant(object, prior, significance)
            (stats < lh$low) | (stats > lh$high)
          })

setMethod('countSignificant', signature('ebCorrelation'),
          function(object, prior, significance, ...) {
            sum(selectSignificant(object, prior, significance))
          })

setMethod('summary', 'ebCorrelation', function(object, prior=1, significance=0.9, ...) {
  lh <- cutoffSignificant(object, prior, significance)
  cat(paste('Call:', as.character(list(object@call)),'\n'))
  cat(paste('Row-by-row correlation analysis with',
            length(object@correlations), 'rows\n\nDistribution of correlations:\n'))
  print(summary(object@correlations))
  cat(paste('With prior =', prior, 'and alpha =', significance, '\n'))
  cat(paste('\tthe upper tail contains', sum(object@correlations > lh$high),
            'values above', lh$high, '\n'))
  cat(paste('\tthe lower tail contains', sum(object@correlations < lh$low),
            'values below', lh$low, '\n'))
})
