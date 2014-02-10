# Copyright (C) Kevin R. Coombes, 2007-2012

setClass('ClusterTest',
         representation(call='call',
                        result='matrix'))

.getUpperTriValues <- function(data) {
  n <- nrow(data)
  joint <- rep(NA, n*(n-1)/2) # room for upper triangular entries
  start <- 1
  for (i in 1:(n-1)) {
    len <- n - i
    joint[start:(start+len-1)] <- data[i, (i+1):n]
    start <- start + len
  }
  joint
}

setMethod('summary', signature(object='ClusterTest'),
          function(object, ...) {
  cat(paste('A', as.character(class(object)), 'object.\n\n'))
  cat(paste('Call:\n\t', as.character(list(object@call))),'\n\n')
  cat('Agreement levels:\n')
  summary(.getUpperTriValues(object@result))
})

setMethod('hist', signature(x='ClusterTest'),
          function(x, ...) {
  agreement <- .getUpperTriValues(x@result)
  hist(agreement, ...)
})

setMethod('image', signature(x='ClusterTest'),
          function(x, dendrogram, ...) {
  if(missing(dendrogram)) {
    h <- heatmap(x@result, symm=TRUE, revC=FALSE, scale='none', ...)
  } else {
    dendrogram <- as.dendrogram(dendrogram)
    h <- heatmap(x@result, Rowv=dendrogram, symm=TRUE, revC=FALSE, scale='none', ...)
  }
  invisible(h)
})

