#########################################################
## NEW GENERIC METHODS

require('methods')
require('oompaBase')

if (!isGeneric("analyze"))
  setGeneric("analyze",
             function(object, ...) { standardGeneric("analyze") }
             )

if (!isGeneric("process"))
  setGeneric("process",
             function(object, action, parameter=NULL)
             { standardGeneric("process") }
             )

if (!isGeneric("plot"))
  setGeneric("plot", function(x, y, ...) standardGeneric("plot"))
if (!isGeneric("hist"))
  setGeneric("hist", function(x, ...) standardGeneric("hist"))
if (!isGeneric("image"))
  setGeneric("image", function(x, ...) standardGeneric("image"))

if (!isGeneric("print"))
  setGeneric("print", function(x, ...) standardGeneric("print"))
if (!isGeneric("summary"))
  setGeneric("summary", function(object, ...) standardGeneric("summary"))

if (!isGeneric("as.data.frame"))
  setGeneric("as.data.frame", function (x,
                                        row.names = NULL,
                                        optional = FALSE, ...) 
             standardGeneric("as.data.frame"))

