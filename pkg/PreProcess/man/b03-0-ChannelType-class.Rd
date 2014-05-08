\name{ChannelType-class}
\alias{ChannelType-class}
\alias{ChannelType}
\alias{print,ChannelType-method}
\alias{show,ChannelType-method}
\alias{summary,ChannelType-method}
\alias{setDesign}
\alias{getDesign}
\docType{class}
\title{Class "ChannelType"}
\description{
  This class represents the "type" of a microarray channel. 
}
\usage{
ChannelType(mk, md, nc, nr, gl, design="")
setDesign(object, design)
getDesign(object)
\S4method{print}{ChannelType}(x, \dots)
\S4method{show}{ChannelType}(object)
\S4method{summary}{ChannelType}(object, \dots)
}
\arguments{
  \item{mk}{character string specifying the name of the manufacturer of the
    microarray (e.g., 'Affymetrix') }
  \item{md}{character string specifying the model of the microarray
    (e.g., 'Hu95A')}
  \item{nc}{scalar integer specifying the number of columns in the array}
  \item{nr}{scalar integer specifying the number of rows in the array}
  \item{gl}{character string specifying the material used to label samples}
  \item{design}{character string containing the name of an object describing
    details about the design of the microarray}
  \item{object}{object of class \code{ChannelType}}
  \item{x}{object of class \code{ChannelType}}
  \item{\dots}{extra arguments for generic or plotting routines}
}
\details{
  Microarrays come in numerous flavors. At present, the two most common
  types are the synthesized oligonucleotide arrays produced by
  Affymetrix and the printed cDNA arrays on glass, which started in Pat
  Brown's lab at Stanford. In earlier days, it was also common to find
  nylon microarrays, with the samples labeled using a radioactive
  isotope. The glass arrays are distinguished from other kinds of
  arrays in that they typically cohybridize two different samples
  simultaneously, using two different fluorescent dyes. The fluorescence
  from each dye is scanned separately, producing two images and thus two
  related sets of data from the same microarray. We refer to these
  parallel data sets within an array as ``channels''.

  An object of the \code{ChannelType} class represents a combination of
  the kind of microarray along with the kind of labeling procedure.
  These objects are intended to be passed around as part of more complex
  objects representing the actual gene expression data collected from
  particular experiments, in order to be able to eventually tie back into
  the description of what spots were laid down when the array was produced.

  The \code{ChannelType} object only contains a high level description
  of the microarray, however. Detailed information about what
  biological material was laid down at each spot on the microarray is
  stored elsewhere, in a ``design'' object. Within a \code{ChannelType}
  object, the design is represented simply by a character string. This
  string should be the name of a separate object containing the detailed
  design information. This implementation allows us to defer the design
  details until later. It also saves space by putting the details in a
  single object instead of copying them into every microarray. Finally,
  it allows that single object to be updated when better biological
  annotations are available, with the benefits spreading immediately to
  all the microarray projects that use that design.
}
\section{Slots}{
  \describe{
    \item{\code{maker}:}{character string specifying the name of the
      manufacturer of the microarray}
    \item{\code{model}:}{character string specifying the model of the
      microarray} 
    \item{\code{nCol}:}{scalar integer specifying number of columns in the array}
    \item{\code{nRow}:}{scalar integer specifying number of rows in the array}
    \item{\code{glow}:}{character string specifying the material used to label
      samples}
    \item{\code{design}:}{character string containing the name of an object
      describing details about the design of the microarray}
  }
}
\section{Methods}{
  \describe{
    \item{print(x, \dots)}{Prints all the information in the object}
    \item{show(object)}{Prints all the information in the object}
    \item{summary(object, \dots)}{Writes out a summary of the object}
  }
}
\value{
  The \code{ChannelType} constructor returns a valid object of the
  class.

  The \code{setDesign} function invisibly returns the \code{ChannelType}
  object on which it was invoked.

  The \code{getDesign} function returns the design object referred to by
  the \code{design} slot in the \code{ChannelType} object. If this string
  does not evaluate to the name of an object, then \code{getDesign} returns
  a \code{NULL} value.
}
\author{
  Kevin R. Coombes \email{krc@silicovore.com},
  P. Roebuck \email{proebuck@mdanderson.org}
}
\seealso{
  \code{\link{Channel}}
}
\examples{
showClass("ChannelType")

x <- ChannelType('Affymetrix', 'oligo', 100, 100, 'fluor')
x
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

