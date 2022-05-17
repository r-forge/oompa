\name{Contribution-class}
\docType{class}
\alias{Contribution}
\alias{Contribution-class}
\alias{[,Contribution,ANY,ANY,ANY-method}%]
\alias{summary,Contribution-method}
\alias{image,Contribution-method}
\alias{heat,Contribution-method}
\alias{getCompositeWeights}
\alias{getAllWeights}
\alias{getTop}
\alias{pickSignificant}

\title{Class \code{"Contribution"}}
\description{
  The \code{Contribution} object class ...
}
\section{Objects from the Class}{
  Objects are defined using the \code{getAllWeights}, \code{getCompositeWeights}, \code{getTop}, or \code{pickSignificant} functions. In the simplest scenario, one would enter an object of class \code{plasma} and any specific parameters associated with the function (see arguments section for more info).
}
\usage{
getCompositeWeights(object, N, M)
getAllWeights(object, N)
getTop(object, N=1)
pickSignificant(object, alpha)
\S4method{summary}{Contribution}(object, M, N, \dots)
\S4method{image}{Contribution,missing}(object, \dots)
\S4method{heat}{Contribution,missing}(object, \dots)
}
\arguments{
  \item{object}{an object of the \code{plasma} class.}
  \item{N}{name of the dataset being modeled.}
  \item{M}{name of the dataset being modeled pairwise with dataset \code{N} in the \code{getCompositeWeights} function.}
  \item{alpha}{level of significance used in the \code{pickSignificant} function.}
  \item{\dots}{other graphical parameters.}
}
\section{Slots}{
  \describe{
      \item{\code{contrib}:}{a matrix of the original variables in dataset \code{N} as rows and the PLS components \code{M} as columns.}
      \item{\code{datasets}:}{a character vector that stores the names of the datasets that were specified for the function.}
  }
}
\section{Methods}{
  \describe{
      \item{\code{\link{summary}}:}{outputs summary statistics for the contributions of dataset \code{N} to components from all datasets in the case of \code{getAllWeights} or dataset \code{M} in the case of \code{getCompositeWeights}.}
      \item{\code{\link{image}}:}{outputs a heatmap of the transposed \code{contrib} matrix.}
      \item{\code{\link{heat}}:}{outputs a clustered heatmap of the \code{contrib} matrix.}
  }
}
\value{
}
\author{Kevin R. Coombes \email{krc@silicovore.com}}
\note{
}
\seealso{}
\examples{
}
\keyword{classes}
