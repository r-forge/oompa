# Copyright (C) Kevin R. Coombes, 2007-2012

##########################################################################
# COMPLETECHANNEL		attributes: name, type, data
#	name is a character string
#	type is an object of class channel.type
#	data is a data.frame containing (at least) vol, bkgd, svol
#
#	print			print all the information
#	summary			summarize the data
#	as.data.frame		fetch the underlying data
#	analyze			estimate density functions
#	plot			plot estimated density functions
#	image			show values 'geographicaly'
#	channelize		what kinds of channels are extracted?

setClass('CompleteChannel',
         representation(name = 'character',
                        type = 'ChannelType',
                        data = 'data.frame',
                        history = 'list'))

# The history object is used to keep a record of the calls that did the
# data processing.
CompleteChannel <- function(name, type, data) {
  cl <- as.character(list(match.call())) # i wish i knew why this works
  new('CompleteChannel', name = name, type = type, data = data,
      history = list(cl))
}

setMethod('print', signature(x='CompleteChannel'),
          function(x, ...) {
	print(x@type)
	print(x@name)
	print(x@data, ...)
})

setMethod('show', signature(object='CompleteChannel'),
          function(object) {
	print(object@type)
	print(object@name)
	print(object@data)
})

setMethod('summary', signature(object='CompleteChannel'),
          function(object, ...) {
	print(object@type)
	print(object@name)
	summary(object@data, ...)
})

setMethod('as.data.frame', signature(x='CompleteChannel'),
          function(x, row.names=NULL, optional=FALSE) {
              x@data
})

# The function analyze.CompleteChannel computes three density functions
# for the given channel: one each for volume, background, and corrected
# volume. These are returned as components of a list.
setMethod('analyze', signature(object='CompleteChannel'),
          function(object, useLog=FALSE, ...) {
  ch <- x@data
  if (useLog) {
    vwid <- 0.4
    bwid <- 0.3
    ch <- data.frame(ch$vol, ch$bkgd, ch$svol)
    ch <- data.frame(logb(ch, 2))
    dimnames(ch)[[2]] <- c('vol', 'bkgd', 'svol')
  } else {
    vwid <- 100
    bwid <- 40
  }
  dvol <- density(ch$vol, n = 300, width = vwid, na.rm = TRUE)
  dbkgd <- density(ch$bkgd, n = 300, width = bwid, na.rm = TRUE)
  dsvol <- density(ch$svol[is.finite(ch$svol)], n = 300,
                   width = vwid, na.rm = TRUE)
  list(dvol = dvol, dbkgd = dbkgd, dsvol = dsvol)
})

# The auxiliary function .f.plot3 simply plots the three density functions 
# produced by analyze.CompleteChannel. It's really not intended for end
# users, but it is also called directly from "complete-slide.R".
.f.plot3 <- function(d, name, ...) {
  plot(d$dvol,  type='l', xlab='Intensity', ylab='Frequency',
       main=paste(name, 'Volume'), ...)
  plot(d$dbkgd, type='l', xlab='Intensity', ylab='Frequency',
       main=paste(name, 'Background'), ...)
  plot(d$dsvol, type='l', xlab='Intensity', ylab='Frequency',
       main=paste(name, 'Corrected Volume'), ...)
}

# The function plot.CompleteChannel combines the computation of channel
# densities with a plotting routine.
setMethod('plot', signature('CompleteChannel', 'missing'),
          function(x, main=x@name, useLog=FALSE, ...) {
            d <- analyze(x, useLog)
            .f.plot3(d, main, ...)
            invisible(x)
})

# When you extract processed channel data from within a complete
# channel, you need to know the correct type for the new object.
# This method is a hook so classes derived from CompleteChannel
# can do the right thing. The return value is used as the first
# argument in a call to "new". Should plan on making this a proper
# method as soon as we get a derived class.

setMethod('channelize', signature(object='CompleteChannel'),
          function(object, ...) { 'Channel' })

setMethod('process', signature('CompleteChannel', 'Processor'),
          function(object, action, parameter=NULL) {
            if (is.null(parameter)) {
              parameter <- action@default
            }
            x <- action@f(object, parameter)
            parent <- deparse(substitute(object))
            result.class <- channelize(object)
            history <- paste(action@description,
                         ' (using object: ', deparse(substitute(action)),
                         ') with parameter = ', parameter, sep='')
            name <- paste(action@name, object@name, sep='')
            new(result.class,
                x=x, parent=parent, type=object@type,
                name=name, history=list(history))
          })

# now we look at various ways to process CompleteChannels.

PROC.BACKGROUND <-
  new('Processor',
      name='background ',
      description='Background intensity',
      default=0,
      f= function(ch, extra=0) { ch@data$bkgd }
      )

PROC.SIGNAL <-
  new('Processor',
      name='signal ',
      description='Foreground intensity',
      default=0,
      f= function(ch, extra=0) { ch@data$vol }
      )


setMethod('image', signature(x='CompleteChannel'),
          function(x, ...) {
  image(process(process(x, PROC.SIGNAL), PROC.LOG.TRANSFORM, 2), ...)
  image(process(process(x, PROC.BACKGROUND), PROC.LOG.TRANSFORM, 2), ...)
  invisible(x)
})

PROC.CORRECTED.SIGNAL <-
  new('Processor',
      name='background-corrected signal ',
      description='Get background-corrected signal',
      default=0,
      f= function(ch, extra=0) { ch@data$svol }
      )

PROC.NEG.CORRECTED.SIGNAL <-
  new('Processor',
      name='background-corrected signal ',
      description='Foreground - Background',
      default=0,
      f= function(ch, extra=0) { ch@data$vol - ch@data$bkgd }
      )

PROC.SD.SIGNAL <-
  new('Processor',
      name='SD ',
      description='Standard Deviation of Foreground',
      default=0,
      f= function(ch, extra=0) { ch@data$SD }
      )

PROC.SIGNAL.TO.NOISE <-
  new('Processor',
      name='signal-to-noise ratio ',
      description='Signal-to-noise ratio',
      default=0,
      f= function(ch, extra=0) { ch@data$SN }
      )

makeDefaultPipeline <- function(ef = PROC.SIGNAL,	ep = 0,
                                nf = PROC.GLOBAL.NORMALIZATION, np = 0,
                                tf = PROC.THRESHOLD, tp = 25,
                                lf = PROC.LOG.TRANSFORM, lp = 2,
                                name='standard pipe',
                                description='my method') {
  base <- ef
  base@default <- ep
  norm <- nf
  norm@default <- np
  thresh <- tf
  thresh@default <- tp
  trans <- lf
  trans@default <- lp
  new('Pipeline',
      proclist=list(base, norm, thresh, trans),
      name=name, description=description)
}


PIPELINE.MDACC.DEFAULT <- makeDefaultPipeline(description='MDACC default')

