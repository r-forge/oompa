\name{MultiplePLSCoxModels-class}
\docType{class}
\alias{MultiplePLSCoxModels}
\alias{MultiplePLSCoxModels-class}
\alias{summary,MultiplePLSCoxModels-method}
\alias{[,MultiplePLSCoxModels,ANY,ANY,ANY-method}%]
\alias{plot,MultiplePLSCoxModels,missing-method}
\alias{validMultipleCoxModels}
\alias{fitCoxModels}
\alias{getSizes}

\title{Class \code{"MultiplePLSCoxModels"}} 
\description{
  The \code{MultiplePLSCoxModels} object class ...
  The \code{validMultipleCoxModels} function checks if each data set contains the same set of samples.
  The \code{fitCoxModels} function fits many plsRcoxmodels and returns an S4 object of class MultiplePLSCoxModels.
  The \code{getSizes} returns a matrix with the list of dataframes of the \code{MultiOmics} object as rownames and columns with NT, cNT, and p-values.
}

\section{Objects from the Class}{

}
\usage{
fitCoxModels(object, timevar, eventvar, eventvalue)
\S4method{summary}{MultiplePLSCoxModels}(object, \dots)
\S4method{plot}{MultiplePLSCoxModels,missing}(x, y, col = c("blue", "red"),
      lwd = 2, xlab = "", ylab = "Fraction Surviving",
       mark.time = TRUE, legloc = "topright", \dots)
}
\arguments{
  \item{object}{an object of class \code{MultiOmics} for fitting the model.}
  \item{timevar}{a column in the \code{MultiOmics} object in the
    \code{outcome} dataframe containing the time-to-event.}
  \item{eventvar}{a column in the \code{MultiOmics} object in the
    \code{outcome} dataframe containing the event.}
  \item{eventvalue}{a character string specifying the value of the event in \code{eventvar}.}
  
  \item{object}{an object of class \code{MultiplePLSCoxModels} for outputting the summary.}
  
  \item{x}{an object of class \code{MultiplePLSCoxModels} for plotting the Kaplan-Meier curves.}
  \item{N}{A character string identifying the data set being modeled.}
  \item{y}{An ignored argrument for the plot method.}
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
      \item{\code{models}:}{A list of \code{SingleModel} objects, one for each assay.}
      \item{\code{timevar}:}{A character matching the name of the column containing the time-to-event.}
      \item{\code{eventvar}:}{A character matching the name of the column containing the event.}
      \item{\code{eventvalue}:}{A character specifying the event in \code{eventvar.}
  }
}
\section{Methods}{
  \describe{
      \item{\code{\link{plot}}:}{Plots Kaplan-Meier curves for each
	omics dataset split into Low Risk and High Risk groups.}
      \item{\code{\link{summary}}:}{Returns a description of the
	\code{MultiplePLSCoxModels} object and the names of the omics
	datasets used to build the model.}
      \item{\code{\link{predict}}:}{spaceholder.}
  }
  }
}
\value{

}
\author{Kevin R. Coombes \email{krc@silicovore.com}}
\note{
}
\seealso{}
\examples{
data("TCGA-ESCA")
MO <- prepareMultiOmics(assemble, Outcome)
splitVec <- rbinom(nrow(Outcome), 1, 0.6)
trainD <- MO[, splitVec == 1]
testD <- MO[, splitVec == 0]

suppressWarnings( firstPass <- fitCoxModels(trainD, timevar = "Days",
                          eventvar = "vital_status", eventvalue = "dead") )
summary(firstPass)
plot(firstPass)
getSizes(firstPass)
}
\keyword{classes}
