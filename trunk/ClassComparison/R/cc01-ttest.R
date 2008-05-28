############################
# vectorized functions to compute row-wise mean and variance quickly
matrixMean <- function(x) {
  x <- as.matrix(x)
  n <- ncol(x)
  x %*% matrix(rep(1, n), nrow=n)/n
}

matrixVar <- function(x, xmean) {
  x <- as.matrix(x)
  n <- ncol(x)
  v <- (((x * x) %*% matrix(rep(1, n), nrow=n)) - n * xmean * xmean)/(n - 1)
  v[v < 0] <- 0
  v
}

matrixT <- function(m,v) {
  v <- v==levels(v)[1]
  am <- matrixMean(m[,v])
  an <- sum(v)
  av <- matrixVar(m[,v], am)
  bm <- matrixMean(m[,!v])
  bn <- sum(!v)
  bv <- matrixVar(m[,!v], bm)
  (am-bm)/sqrt(((an-1)*av + (bn-1)*bv)/(an+bn-2))/sqrt(1/an+1/bn)
}

setClass('MultiTtest',
         representation(t.statistics='numeric',
                        p.values='numeric',
                        groups='character',
                        call='call'))

MultiTtest <- function(data, classes) {
  call <- match.call()
  if(is.logical(classes)) classes <- factor(classes)
  if(inherits(data, 'ExpressionSet')) {
    if(is.character(classes)) {
      classes <- as.factor(pData(data)[,classes])
    }
    data <- exprs(data)
  }
  t.statistics <- as.vector(matrixT(data, classes))
  p.values <- sapply(t.statistics, function(tv, df){
    2*(1-pt(abs(tv), df))
  }, length(classes)-2)
  new('MultiTtest',
      t.statistics=t.statistics,
      p.values=p.values,
      groups=levels(classes),
      call=call)
}

setMethod('summary', 'MultiTtest', function(object,...) {
  cat(paste('Row-by-row two-sample t-tests with',
            length(object@t.statistics), 'rows\n'))
  cat(paste('Positive sign indicates an increase in class:',
            object@groups[1],'\n\n'))
  cat(paste('Call:', as.character(list(object@call)),'\n\nT-statistics:\n'))
  print(summary(object@t.statistics))
  cat('\nP-values:\n')
  summary(object@p.values)
})

setMethod('hist', 'MultiTtest', function(x, xlab='T Statistics',
                                         main=NULL, ...) {
  hist(x@t.statistics, xlab=xlab, main=main, ...)
})

setMethod('plot', signature('MultiTtest', 'missing'), function(x, y,
                  ylab='T Statistics', ...) {
  plot(x@t.statistics, ylab=ylab, ...)
})

setMethod('plot', signature('MultiTtest', 'ANY'), function(x, y,
                  xlab='T Statistics',
                  ylab=deparse(substitute(y)), ...) {
  plot(x@t.statistics, y, xlab=xlab, ylab=ylab, ...)
})

