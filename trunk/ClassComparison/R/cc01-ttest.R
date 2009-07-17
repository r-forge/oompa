############################
# vectorized functions to compute row-wise mean and variance quickly
matrixMean <- function(x, na.rm=FALSE) {
  x <- as.matrix(x)
  Nper <- n <- ncol(x)
  if (na.rm) {
    Nper <- apply(x, 1, function(y) sum(!is.na(y)))
    x[is.na(x)] <- 0
  }
  x %*% matrix(rep(1, n), nrow=n)/Nper
}

matrixVar <- function(x, xmean, na.rm=FALSE) {
  x <- as.matrix(x)
  Nper <- n <- ncol(x)
  if (na.rm) {
    Nper <- apply(x, 1, function(y) sum(!is.na(y)))
    x[is.na(x)] <- 0
  }
  v <- (((x * x) %*% matrix(rep(1, n), nrow=n)) - Nper * xmean * xmean)/(Nper - 1)
  v[v < 0] <- 0
  v
}

matrixT <- function(m, v, na.rm=FALSE) {
  v <- v==levels(v)[1]
  if (na.rm) {
    an <- apply(m[,v], 1, function(y) sum(!is.na(y)))
    bn <- apply(m[,!v], 1, function(y) sum(!is.na(y)))
  } else {
    an <- sum(v)
    bn <- sum(!v)
  }
  am <- matrixMean(m[,v],     na.rm=na.rm)
  av <- matrixVar(m[,v],  am, na.rm=na.rm)
  bm <- matrixMean(m[,!v],    na.rm=na.rm)
  bv <- matrixVar(m[,!v], bm, na.rm=na.rm)
  (am-bm)/sqrt(((an-1)*av + (bn-1)*bv)/(an+bn-2))/sqrt(1/an+1/bn)
}

setClass('MultiTtest',
         representation(t.statistics='numeric',
                        p.values='numeric',
                        df='numeric',
                        groups='character',
                        call='call'))

MultiTtest <- function(data, classes, na.rm=TRUE) {
  call <- match.call()
  if(is.logical(classes)) classes <- factor(classes)
  if(inherits(data, 'ExpressionSet')) {
    if(is.character(classes)) {
      classes <- as.factor(pData(data)[,classes])
    }
    data <- exprs(data)
  }
  t.statistics <- as.vector(matrixT(data, classes, na.rm=na.rm))
  df <- apply(data, 1, function(x) sum(!is.na(x))-2)
#  p.values <- sapply(t.statistics, function(tv, df){
#    2*(1-pt(abs(tv), df))
#  }, length(classes)-2)
  p.values <- 2*(1-pt(abs(t.statistics), df))
  new('MultiTtest',
      t.statistics=t.statistics,
      p.values=p.values,
      df=df,
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

setMethod("as.data.frame", "MultiTtest", function (x, row.names = NULL,
                                                   optional = FALSE, ...) {
  data.frame(Tstats=x@t.statistics,
             Pvalues=x@p.values,
             DF=x@df)
})
