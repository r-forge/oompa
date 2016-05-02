# Copyright (C) Kevin R. Coombes, 2007-2016

##
## S3 or non-generics converted to S4
##

if (!isGeneric("as.data.frame"))
  setGeneric("as.data.frame",
             function(x, row.names=NULL, optional=FALSE, ...) 
                 standardGeneric("as.data.frame"))

if (!isGeneric("hist"))
  setGeneric("hist",
             function(x, ...) standardGeneric("hist"))

if (!isGeneric("image"))
  setGeneric("image",
             function(x, ...) standardGeneric("image"))

if (!isGeneric("plot"))
  setGeneric("plot",
             function(x, y, ...) standardGeneric("plot"))

if (!isGeneric("print"))
  setGeneric("print",
             function(x, ...) standardGeneric("print"))

if (!isGeneric("summary"))
  setGeneric("summary",
             function(object, ...) standardGeneric("summary"))

if (!isGeneric("predict"))
  setGeneric("predict", function(object, ...) standardGeneric("predict"))

if (!isGeneric("screeplot"))
  setGeneric("screeplot", function(x, ...)
             standardGeneric("screeplot"))

if (!isGeneric("pltree"))
  setGeneric("pltree", function(x, ...) standardGeneric("pltree"))
