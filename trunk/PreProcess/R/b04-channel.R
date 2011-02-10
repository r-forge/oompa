##########################################################################
# CHANNELS		attributes: parent, name, type, x
#
#	channel		constructor
#	plot		values as a function of position
#	hist		values
#	summary		underlying vector
#	print		underlying vector
#	image		geographic display of values

setClass('Channel',
         representation(parent = 'character',
                        name = 'character',
                        type = 'ChannelType',
                        x = 'numeric',
                        history = 'list'))

# The history object is used to keep a record of the calls that did the
# data processing.
Channel <- function(parent, name, type, vec) {
  cl <- as.character(list(match.call())) # i wish i knew why this works
  new('Channel', parent=parent, name=name, type=type, x=vec,
      history=list(cl))
}

# The plot method for a channel produces a graph of the values against
# the position.
setMethod('plot', signature('Channel', 'missing'), function(x, ...) {
  plot(x@x, xlab='Position', ylab=x@name, main=x@parent, ...)
  invisible(x)
})

setMethod('hist', 'Channel', function(x, breaks=67,
                                      xlab=x@name,
                                      main=x@parent, ...) {
  hist(x@x, breaks=breaks, xlab=xlab, main=main, ...)
  invisible(x)
})

# The summary method for a channel summarizes the underlying vector.
setMethod('summary', 'Channel', function(object, ...) {
  cat(paste(object@name, ', a microarray channel object\n', sep=''))
  parent <- object@parent
  if (parent == '') parent <- 'NA'
  cat(paste('Parent object:', parent, '\n'))
  summary(object@type)
  cat('History:\n')
  lapply(object@history, function(x) {
    cat(paste('\t', x, '\n'))
  })
  cat('\n')
  summary(object@x, ...)
})

# The print method for a channel prints the underlying vector.
setMethod('print', 'Channel', function(x, ...) {
  print(x@type)
  print(x@x, ...)
})

# The image method for a slide produces a cartoon of the geographically
# arranged values for each channel.
# Quantiles are reported in the subtitle of each figure.
setMethod('image', 'Channel', function(x, main=x@name, sub=NULL, ...) {
  my.show <- function(x) {
    qq <- quantile(x, c(0.25, 0.5, 0.75), na.rm=TRUE)
    paste('median = ', format(qq[2], digits=4),
          ', IQR = ', format(qq[3]-qq[1], digits=4))
  }
  my.cartoon <- function(v, ncol, ...) {
    if (!is.matrix(v)) {
      v <- t(matrix(v, ncol=ncol, byrow=TRUE))
    }
    v <- fliplr(v)
    image(v, ...)
  }
  nc <- x@type@nCol
  nr <- x@type@nRow
  if (is.null(sub)) {
    sub <- my.show(x@x)
  }
  my.cartoon(x@x[1:(nc*nr)], ncol=nc, main=main, sub=sub)
  invisible(x)
}) 


