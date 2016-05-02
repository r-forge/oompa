# Copyright (C) Kevin R. Coombes, 2007-2016

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

# vectorized functions to get row-by-row t-tests (with variants) quickly

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

## paired t-test
matrixPairedT <- function (m, v, pf)
{
  v <- v == levels(v)[1]
  m <- m[,order(v,pf)] # fix according to Hui Yao
  am <- m[,v]
  bm <- m[,!v]
  pd <- am-bm                   ## paired difference
  pdm <- matrixMean(pd)         ## mean and...
  pdv <- matrixVar(pd,pdm)      ## variance of paired diffs
  n <- ncol(am)
  pdm/(sqrt(pdv/(n)))
}

## unequal variance t-test
matrixUnequalT <- function (m, v) {
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


