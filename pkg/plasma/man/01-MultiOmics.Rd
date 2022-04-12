\name{MultiOmics-class}
\docType{class}
\alias{MultiOmics}
\alias{MultiOmics-class}
\alias{summary,MultiOmics-method}
\alias{[,MultiOmics,ANY,ANY,ANY-method}%]
\alias{plot,MultiOmics,missing-method}
\alias{validMultiOmics}
\alias{prepareMultiOmics}

\title{Class \code{"MultiOmics"}}
\description{
  The \code{MultiOmics} object class ...
  The \code{prepareMultiOmics} function returns a new object of \code{MultiOmics} class for use in \code{fitCoxModel}.
}
\section{Objects from the Class}{
  Objects should be defined using the \code{prepareMultiOmics} constructor. In
  the simplest case, you enter two objects: a list of dataframes and a dataframe of clinical outcomes.
}
\usage{
prepareMultiOmics(datalist, outcome)
\S4method{summary}{MultiOmics}(object, \dots)
\S4method{plot}{MultiOmics,missing}(x, y, \dots)
}
\arguments{
  \item{datalist}{a list of dataframes formatted to have sample names as column indexes and variable names as row indexes}
  \item{outcome}{a dataframe of clinical outcomes formatted to have
    sample names as row indexes and variable names as column indexes}
  \item{object}{An object of the \code{MultiOmics} class.}
  \item{x}{An object of the \code{MultiOmics} class.}
  \item{y}{Nothing; ignored.}
  \item{\dots}{Extra graphical or other parameters.}
}
\section{Slots}{
  \describe{
      \item{\code{data}:}{A list of dataframes with variables as rows and samples as columns of uniform length.}
      \item{\code{outcome}:}{A dataframe of clinical outcomes with variables as columns and samples as rows.}
  }
}
\section{Methods}{
  \describe{
      \item{\code{\link{plot}}:}{Produces a visual representation of the  dimensionalities of each dataframe in datalist. Gray areas correspond to missing samples.}
      \item{\code{\link{summary}}:}{Produces summary tables corresponding to datasets and outcomes.}
  }
}
\value{
The \code{prepareMultiOmics} function returns a new object of the \code{MultiOmics} class.
}
\author{Kevin R. Coombes \email{krc@silicovore.com}}
\note{
}
\seealso{}
\examples{
data("TCGA-ESCA")
prep_res <- prepareMultiOmics(datalist = assemble, outcome = Outcome)

plot(prep_res)
summary(prep_res)
}
\keyword{classes}
