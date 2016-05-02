# Copyright (C) Kevin R. Coombes, 2007-2016

################################################################
# COLOR.CODING		attributes: v, color, mark
#
#    ColorCoding	constructor from a logical vector and two scalars

setClass('ColorCoding',
         representation(v='logical', color='vector', mark='vector'),
         prototype=list(v=T, color='red', mark=16))

ColorCoding <- function(v, color, mark=1) {
  val <- new('ColorCoding', v=v, color=color, mark=mark)
}

################################################################
# COLOR.CODED.PAIR		attributes: x, y, ccl
#
#	ColorCodedPair	constructor from two vectors and a list
#				that contains ColorCodings
#	plot			scatter plot with colors

setClass('ColorCodedPair',
         representation(x='numeric', y='numeric', colorCodingList='list'))

ColorCodedPair <- function(x, y, ccl) {
  val <- new('ColorCodedPair', x=x, y=y, colorCodingList=ccl)
}

setMethod('plot', signature(x='ColorCodedPair', y='missing'),
          function(x, ...) {
  ob <- x
  myargs <- list(...)
  cex <- myargs$cex
  if(is.null(cex)) cex <- par('cex')
  plot(ob@x, ob@y, type='n', ...)
  if (is.list(ob@colorCodingList)) {
    lapply(ob@colorCodingList, function(cc, x, y) {
      if (sum(cc@v) > 0) {
        points(x[cc@v], y[cc@v], col=cc@color, pch=cc@mark, cex=cex)
      }
    }, ob@x, ob@y)
  } else if (is(ob@colorCodingList, 'ColorCoding')) {
    v <- ob@colorCodingList@v
    if (sum(v) > 0) {
      points(ob@x[v], ob@y[v], col=ob@colorCodingList@color,
             pch=ob@colorCodingList@mark, cex=cex)
    }
  }
  invisible(ob)
}
)

colorCode <- function(fac, colorScheme=1:length(levels(fac)), mark=1) {
  lapply(1:length(levels(fac)), function(i, fac, cs, m) {
    lev <- levels(fac)[i]
    ColorCoding(fac==lev, cs[i], m)
  }, fac, colorScheme, mark)
}
