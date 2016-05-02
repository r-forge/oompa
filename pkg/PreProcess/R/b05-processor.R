# Copyright (C) Kevin R. Coombes, 2007-2016

# A "processor" represents a function that acts on the data of a
# some object to process it in some way. The result is always
# another related object, which should record some history about
# exactly how it was processed.

setClass('Processor',
         representation(f = 'function',
                        default = 'numeric or NULL',
                        name = 'character',
                        description = 'character'))

# To demonstrate, we first show how to process channel objects. Each
# channel basically represents a single vector of data, together with
# enough meta information to understand where it came from.

setMethod('process', signature('Channel', 'Processor'),
          function(object, action, parameter = NULL) {
            if (is.null(parameter)) {
              parameter <- action@default
            }
            object@x <- action@f(object@x, parameter)
            what <- paste(action@description,
                          ' (using object: ', deparse(substitute(action)),
                          ') with parameter = ', parameter, sep = '')
            object@history <- append(object@history, what)
            object@name <- paste(action@name, object@name, sep = '')
            object
          })

setMethod('summary', signature(object='Processor'),
          function(object, ...) {
  cat(paste('A processor object ', object@name, '\nAction: ',
            object@description, '\nDefault parameter: ',
            object@default, '\n',
            sep = ''))
})

PROC.SUBTRACTOR <-
  new('Processor',
      name = '',
      description = 'Subtracted a global constant',
      default = 0,
      f = function(x, v) { x - v}
      )

PROC.THRESHOLD <-
  new('Processor',
      name = '',
      description = 'Truncated below',
      default = 0,
      f = function(x, th) {
        x[x < th] <- th;
        x
      })

PROC.GLOBAL.NORMALIZATION <-
  new('Processor',
      name = 'normalized ',
      description = 'Global normalization',
      default = 0,
      f = function(x, n) {
        if (!(n>0)) {
          n <- quantile(x, 0.75, na.rm=TRUE)/1000
        }
        x/n
      })

PROC.LOG.TRANSFORM <-
  new('Processor',
      name = 'log ',
      description = 'Log transformation',
      default = 2,
      f = logb)

PROC.MEDIAN.EXPRESSED.NORMALIZATION <-
  new('Processor',
      name = 'normalized ',
      description = 'Normalized to the median of expressed values',
      default = 0,
      f = function(x, n) {
        if(!(n > 0)) {
          n <- median(x[x > 0])/1000
        }
        x/n
      })

PROC.SUBSET.NORMALIZATION <-
  new('Processor',
      name = 'normalized ',
      description = 'normalized to the median of a subset',
      default = 0,	# means use the global median
      f = function(x, n) {
        if (length(n) < 2) {
          n <- rep(T, length(x))
        }	
        x/median(x[n])*1000
      })

PROC.SUBSET.MEAN.NORMALIZATION <-
  new('Processor',
      name = 'normalized ',
      description = 'normalized to the mean of a subset',
      default = 0,	# means use the global mean
      f = function(x, n) {
        if (length(n) < 2) {
          n <- rep(T, length(x))
        }	
        x/mean(x[n])*1000
      })


