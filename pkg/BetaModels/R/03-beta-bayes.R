# use empirical bayes to determine significance of correlation coefficients

setClass("ebCorrelation",
         slots = c(nObservations='numeric'),
         contains = "MultiWilcoxonTest")

# input:
#  ss = is a vector of computed correlation coefficients
# nObs = number of observations used for each cor coef
#  nPoints = number of points at which to fit the distribution
ebCorrelation <- function(ss, nObs, nPoints=500) {
  call <- match.call()
  # start by estimating the empirical distribution
  eps <- 1/(2*nPoints)
  xvals <- seq(-1-eps, 1+eps, length=nPoints)
  pdf <- hist(ss, breaks=xvals, plot = FALSE)$density
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
  xclip <- xvals[xvals >= min(X) & xvals <= max(X)]
  YP <- predict(Z, data.frame(X = xclip))
  YP <- c(rep(NA, sum(xvals < min(X))),
          YP, rep(NA, sum(xvals > max(X))))
  # convert back to the original scale
  unravel <- exp(YP+log(theo))
  unravel[is.na(unravel)] <- min(unravel, na.rm = TRUE)/10
  new("ebCorrelation",
      nObservations=nObs, # new slots
      xvals=xvals, statistics=ss, pdf=pdf, theoretical.pdf=theo,
      unravel=unravel, call=call,  # old slots
      groups = character(0)) # useles old slots
}


setMethod('hist', 'ebCorrelation', function(x,
           xlab='Correlation',
           ylab='Prob(Different | Y)', main='', ...) {
  callNextMethod(x, xlab = xlab, ylab = ylab, main = main, ...)
  invisible(x)
})

setMethod('plot', signature('ebCorrelation', 'missing'), function(x,
                  prior=1, significance=0.9, ylim=c(-0.5, 1),
                  xlab='Correlation',
                  ylab='Prob(Unusual | Rho)', ...) {
  callNextMethod(x, prior = prior, significance = significance,
                 ylim = ylim, xlab=xlab, ylab=ylab, ...)
  invisible(x)
})

setMethod('summary', 'ebCorrelation', function(object, prior=1, significance=0.9, ...) {
  lh <- cutoffSignificant(object, prior, significance)
  cat(paste('Call:', as.character(list(object@call)),'\n'))
  cat(paste('Row-by-row correlation analysis with',
            length(object@statistics), 'rows\n\nDistribution of correlations:\n'))
  print(summary(object@statistics))
  cat(paste('With prior =', prior, 'and alpha =', significance, '\n'))
  cat(paste('\tthe upper tail contains', sum(object@statistics > lh$high),
              'values above', lh$high, '\n'))
  cat(paste('\tthe lower tail contains', sum(object@statistics < lh$low),
            'values below', lh$low, '\n'))
})
