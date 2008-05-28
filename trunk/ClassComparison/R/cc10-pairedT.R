## paired t-tests, vectorized
## basic code for the paired t-0test was contributed by Roland Bassett.

## utility function
matrixpairedT <- function (m, v, pf)
{
  v <- v == levels(v)[1]
  m <- m[,order(pf,v)]
  am <- m[,v]
  bm <- m[,!v]
  pd <- am-bm                   ## paired difference
  pdm <- matrixMean(pd)         ## mean and...
  pdv <- matrixVar(pd,pdm)      ## variance of paired diffs
  n <- ncol(am)
  pdm/(sqrt(pdv/(n)))
}

## class definition
setClass('MultiTtestPaired',
         representation('MultiTtest'))

MultiTtestPaired <- function (data, classes, pairing)
{
  call <- match.call()
  if (is.logical(classes))
    classes <- factor(classes)
  if (inherits(data, "ExpressionSet")) {
    if (is.character(classes)) {
      classes <- as.factor(pData(data)[, classes])
    }
    data <- exprs(data)
  }
  t.statistics <- as.vector(matrixpairedT(data, classes, pairing))
  p.values <- sapply(t.statistics, function(tv, df) {
    2 * (1 - pt(abs(tv), df))
  }, length(unique(pairing)) - 1)
  new("MultiTtestPaired", t.statistics = t.statistics, p.values = p.values,
      groups = levels(classes), call = call)
}

setMethod('summary', 'MultiTtestPaired', function(object,...) {
  cat('Results of a paired t-test\n')
  callNextMethod(...)
})

##########################################
## unequal variance t tests

## utility function
matrixunequalT <- function (m, v) {
  v <- v==levels(v)[1]
  am <- matrixMean(m[,v])
  an <- sum(v)
  av <- matrixVar(m[,v], am)
  bm <- matrixMean(m[,!v])
  bn <- sum(!v)
  bv <- matrixVar(m[,!v], bm)
  tt <- (am-bm)/sqrt(av/an + bv/bn)
  u <- bv/av
  df <- trunc((1/an + u/bn)^2/(1/(an^2*(an-1)) + u^2/(bn^2*(bn-1))))
  list(tt=tt, df=df)
}


## class definition
setClass('MultiTtestUnequal',
         representation('MultiTtest',
                        df='numeric'))

MultiTtestUnequal <- function (data, classes)
{
  call <- match.call()
  if (is.logical(classes))
    classes <- factor(classes)
  if (inherits(data, "ExpressionSet")) {
    if (is.character(classes)) {
      classes <- as.factor(pData(data)[, classes])
    }
    data <- exprs(data)
  }
  temp <- matrixunequalT(data,classes)
  t.statistics <- as.vector(temp$tt)
  df <- as.vector(temp$df)
  p.values <- 2 * (1 - pt(abs(t.statistics), df))
  new("MultiTtestUnequal", t.statistics = t.statistics, p.values = p.values,
      groups = levels(classes), call = call, df=df)
}

setMethod('summary', 'MultiTtestUnequal', function(object,...) {
  cat('Results of an unequal variance t-test\n')
  callNextMethod(...)
})

