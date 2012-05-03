\name{ChannelType-class}
\alias{ChannelType-class}
\alias{ChannelType}
\alias{print,ChannelType-method}
\alias{summary,ChannelType-method}
\alias{setDesign}
\alias{getDesign}

\title{Class "ChannelType"}
\description{
  This class represents the "type" of a microarray channel. 
}
\usage{
ChannelType(mk, md, nc, nr, gl, design = "")
setDesign(object, design)
getDesign(object)
\S4method{print}{ChannelType}(x,...)
\S4method{summary}{ChannelType}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mk}{ A string; the name of the manufacturer of the microarray
    (e.g., 'Affymetrix') }
  \item{md}{ A string; the model of the microarray (e.g., 'Hu95A')}
  \item{nc}{ An integer; the number of columns in the array}
  \item{nr}{ An integer; the number of rows in the array}
  \item{gl}{ A string; the material used to label samples}
  \item{design}{ A string containing the name of an object describing
    details about the design of the microarray}
  \item{object}{A \code{ChannelType} object}
  \item{x}{A \code{ChannelType} object}
  \item{\dots}{Additional arguments are as in the underlying generic methods.}
}
\details{
  Microarrays come in numerous flavors. At present, the two most common
  types are the synthesized oligonucleotide arrays produced by
  Affymetrix and the printed cDNA arrays on glass, which started in Pat
  Brown's lab at Stanford. In earlier days, it was also common to find
  nylon microarrays, with the samples labeled using a radioactive
  isotope.  The glass arrays are distinguished from other kinds of
  arrays in that they typically cohybridize two different samples
  simultaneously, using two different fluorescent dyes.  The fluorescence
  from each dye is scanned separately, producing two images and thus two
  related sets of data from the same microarray. We refer to these
  parallel data sets within an array as ``channels''.

  An object of the \code{ChannelType} class represents a combination of
  the kind of microarray along with the kind of labeling
  procedure.  These objects are intended to be passed around as part of
  more complex objects representing the actual gene expression data
  collected from particular experiments, in order to be able to
  eventually tie back into the description of what spots were laid down
  when the array was produced.

  The \code{ChannelType} object only contains a high level description
  of the microarray, however. Detailed information about what
  biological material was laid down at each spot on the microarray is
  stored elsewhere, in a ``design'' object. Within a \code{ChannelType}
  object, the design is represented simply by a character string. This
  string should be the name of a separate object containing the detailed
  design information. This implementation allows us to defer the design
  details until later.  It also saves space by putting the details in a
  single object instead of copying them into every microarray. Finally,
  it allows that single object to be updated when better biological
  annotations are available, with the benefits spreading immediately to
  all the microarray projects that use that design.
}
\section{Slots}{
  \describe{
    \item{\code{maker}:}{ A string; the name of the manufacturer of the
      microarray (e.g., 'Affymetrix') }
    \item{\code{model}:}{ A string; the model of the microarray (e.g.,
      'Hu95A')} 
    \item{\code{nCol}:}{ An integer; the number of columns in the array}
    \item{\code{nRow}:}{ An integer; the number of rows in the array}
    \item{\code{glow}:}{ A string; the material used to label samples}
    \item{\code{design}:}{ A string containing the name of an object
      describing details about the design of the microarray}
  }
}
\section{Methods}{
  \describe{
    \item{print(x, ...)}{Prints all the information in the object}
    \item{summary(object, ...)}{Writes out a summary of the object}
  }
}
\value{
  The \code{ChannelType} constructor returns a valid object of the
  class.

  The \code{setDesign} function invisibly returns the \code{ChannelType}
  object on which it was invoked.

  The \code{getDesign} function returns the design object referred to by
  the \code{design} slot in the \code{ChannelType} object.  If this
  string does not evaluate to the name of an object, then
  \code{getDesign} returns a \code{NULL} value.
}
\author{
  Kevin R. Coombes \email{kcoombes@mdanderson.org}
}
\seealso{
  \code{\link{Channel}}
}
\examples{
x <- ChannelType('Affymetrix', 'oligo', 100, 100, 'fluor')
print(x)

summary(x)

y <- setDesign(x, 'fake.design')
print(y)
summary(y)
d <- getDesign(y)
d

rm(d, x, y) # cleanup
}
\keyword{classes}
\keyword{manip}

