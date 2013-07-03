# Copyright (C) Kevin R. Coombes, 2007-2012

###
### ALLGENERICS.R
###


##
## Brand new generics
##

if (!isGeneric("analyze"))
  setGeneric("analyze",
             function(object, ...) standardGeneric("analyze"))

if (!isGeneric("channelize"))
  setGeneric("channelize",
             function(object, ...) standardGeneric("channelize"))

if (!isGeneric("process"))
  setGeneric("process",
             function(object, action, parameter=NULL)
                 standardGeneric("process"))

