# Copyright (C) Kevin R. Coombes, 2007-2016

#########################################################
## NEW CLASS UNION

## This is a kludge. You cannot define this class union in
## a package that also defines objects that use the class union.
## If you try to load the resulting library into an empty R
## workspace, it fails. According to a message on R-devel from
## John Chambers on 15 October 2004, "The fix is simple, though
## I may not get a chance to test & commit it for a couple of days."
## I have no idea whether it was actually fixed, since there was
## no obvious follow-up message. I know that the problem occured
## trying to load the library into R 1.9.0 in June 2005 (on a machine
## that I hadn't upgraded).

## So the fix is to make a separate package that defines the class union,
## and require that package in the one that wants to use it.

# The 'numeric or NULL' class is used for the 'default' slot in
# the 'Processor'' class. This is needed to allow either numeric
# vectors or NULL values to be supplied and stored in the objects.

require(methods)

setClassUnion('numeric or NULL', c('numeric', 'NULL'))
