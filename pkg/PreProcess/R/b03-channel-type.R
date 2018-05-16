# Copyright (C) Kevin R. Coombes, 2007-2016

#################################################################
# CHANNEL.TYPE		attributes: make, model, ncol, nrow, glow
#

setClass('ChannelType',
         slots = c(maker='character',
                   model='character',
                   nCol='numeric',
                   nRow='numeric',
                   glow='character',
                   design='character'))

ChannelType <- function(mk, md, nc, nr, gl, design='') {
  new('ChannelType',  maker=mk, model=md,
      nCol=nc, nRow=nr, glow=gl, design=design)
}

setMethod('print', signature(x='ChannelType'),
          function(x, ...) {
  cat(paste('Microarray type:', x@maker, x@model, '\n'))
  cat(paste('Labeled with:', x@glow), '\n')
})

setMethod('show', signature(object='ChannelType'),
          function(object) {
  cat(paste('Microarray type:', object@maker, object@model, '\n'))
  cat(paste('Labeled with:', object@glow), '\n')
})

setMethod('summary', signature(object='ChannelType'),
          function(object, ...) {
  show(object)
  cat("Design size:", object@nCol, "by", object@nRow, "\n")
  cat("Design information object:", object@design, "\n")
})

# A channel type should know how to get the list of genes that
# were printed in the different locations on the microarray.
# Naturally, that's more complicated than the simple geometry.
# We handle this problem by using the name of a design object
# as a poor man's reference. Of course, a design object will
# probably just be a data.frame with the gene locations as the
# row names, but we can defer that decision for a little while.

getDesign <- function(object) {
  des <- object@design
  x <- try(eval(as.name(des)), silent=TRUE)
  if (is(x, 'try-error')) {
    warning('object does not contain a valid design');
    x <- NULL
  }
  x
}

setDesign <- function(object, design) {
  if (!is(design, 'character')) {
    design <- deparse(substitute(design))
  }
  object@design <- design
  invisible(object)
}
