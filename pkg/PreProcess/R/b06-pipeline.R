# Copyright (C) Kevin R. Coombes, 2007-2016

# A "pipeline" represents a series of processors that should be applied
# in order. You should think of a pipeline as a completely defined set
# of transformations that is applied uniformly to every microarray in
# the set.

setClass('Pipeline',
         slots = c(proclist = 'list',
                   name = 'character',
                   description = 'character'))

setMethod('process', signature('ANY', 'Pipeline'),
          function(object, action, parameter = NULL) {
            # ignore the parameter; pipelines are parameter-free
            what <- paste(action@description,
                          ' (using pipeline: ',
                          deparse(substitute(action)), ')', sep = '')
            object@history <- append(object@history, what)
            for (proc in action@proclist) {
              object <- process(object, proc)
            }
            object
          })

setMethod('summary', signature(object='Pipeline'),
          function(object, ...) {
  cat(paste('A pipeline object:', object@name, '\nGoal:',
            object@description, '\nConcatenates',
            length(object@proclist), 'steps:\n'))
  lapply(object@proclist, function(x) {cat('\n'); summary(x)})
  invisible(object)
})

PIPELINE.STANDARD <- new('Pipeline',
                         proclist = list(PROC.GLOBAL.NORMALIZATION,
                                         PROC.THRESHOLD,
                                         PROC.LOG.TRANSFORM),
                         name = 'standard processor',
                         description = 'Default channel processing')
