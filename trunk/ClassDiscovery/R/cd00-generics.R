require("stats")
require("cluster")

if (!isGeneric("predict"))
  setGeneric("predict", function(object, ...) standardGeneric("predict"))

if (!isGeneric("pltree"))
  setGeneric("pltree", function(x, ...) standardGeneric("pltree"))

if (!isGeneric("screeplot"))
  setGeneric("screeplot", function(x, ...)
             standardGeneric("screeplot"))


