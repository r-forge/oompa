\name{plasmaPredictions-class}
\docType{class}
\alias{plasmaPredictions}
\alias{plasmaPredictions-class}
\alias{plot,plasmaPredictions,missing-method}

\title{Class \code{"plasmaPredictions"}}
\description{
  The \code{plasmaPredictions} object class is returned when running the \code{predict} method on an object of class \code{plasma}.
}
\section{Objects from the Class}{
}
\usage{
\S4method{plot}{plasmaPredictions,missing}(x, y, \dots)
}
\arguments{
  \item{x}{An object of the \code{plasmaPredictions} class for plotting the Kaplan-Meier curves.}
  \item{y}{An ignored argument for the plot method.}
  \item{col}{A vector of color specifications. Default is c(\dQuote{blue}, \dQuote{red}).}
  \item{lwd}{A vactor specifying the line width. Default is \dQuote{2}.}
  \item{xlab}{A character string to label the x-axis. Default is \dQuote{}.}
  \item{ylab}{A character string to label the y-axis. Default is \dQuote{Fraction Surviving}.}
  \item{mark.time}{A logical value; should tickmarks indicate censored
    data? Default is TRUE.}
  \item{legloc}{A character string indicating where to put the legend. Default is \dQuote{topright}.}
  \item{\dots}{Other graphical parameters.}
}
\section{Slots}{
  \describe{
      \item{\code{meanPredictions}:}{A matrix with samples as rows and factors as columns.}
      \item{\code{riskDF}:}{A data.frame of the test data.}
      \item{\code{riskModel}:}{A coxph object...}
      \item{\code{splitModel}:}{A coxph object...}
      \item{\code{SF}:}{A survfit object...}      
  }
}
\section{Methods}{
  \describe{
      \item{\code{\link{plot}}:}{Produces Kaplan-Meier curves for the low risk and high risk groups.}
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
