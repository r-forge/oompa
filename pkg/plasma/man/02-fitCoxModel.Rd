\name{SingleModel-class}
\docType{class}
\alias{fitSingleModel}
\alias{SingleModel-class}
\alias{predict,SingleModel-method}
\alias{summary,SingleModel-method}
\alias{plot,SingleModel,missing-method}
\alias{fitSingleModel}

\title{Class \code{"SingleModel"}}
\description{
  The \code{SingleModel} object class ...
  The \code{fitSingleModel} function takes in an object of
  \code{MultiOmics} class and returns a new object of
  \code{MultiplePLSCoxModels} class. 
}
\section{Objects from the Class}{

}
\usage{
fitSingleModel(object, N, timevar, eventvar, eventvalue)
\S4method{summary}{SingleModel}(object, \dots)
\S4method{plot}{SingleModel,missing}(x, y, col = c("blue", "red"),
      lwd = 2, xlab = "", ylab = "Fraction Surviving",
       mark.time = TRUE, legloc = "topright", \dots)
}
\arguments{
  \item{object}{an object of class \code{MultiOmics} for fitting the model.}
  \item{N}{A character string identifying the data set being modeled.}
  \item{timevar}{a column in the \code{MultiOmics} object in the
    \code{outcome} dataframe containing the time-to-event.}
  \item{eventvar}{a column in the \code{MultiOmics} object in the
    \code{outcome} dataframe containing the event.}
  \item{eventvalue}{a character string specifying the value of the event.}
  
  \item{x}{an object of class \code{plsRcoxmodel} for plotting the Kaplan-Meier curves.}
  \item{y}{An ignored argrument for the plot method.}
  \item{col}{A vector of color specifications.}
  \item{lwd}{A vactor specifying the line width.}
  \item{xlab}{A character string to label the x-axis.}
  \item{ylab}{A character string to label the y-axis.}
  \item{mark.time}{A logical value; should tickmarks indicate censored
    data?}
  \item{legloc}{A character string indicating where to put the legend.}
  \item{\dots}{Other graphical parameters.}
}
\section{Slots}{
  \describe{
      \item{\code{plsmod}:}{Object of class \code{plsRcoxmodel}
	containing the fitted model.}
      \item{\code{Xout}:}{Object of type \code{data.frame} containing
	the original \code{outcome} dataframe and additional columns for
	"Risk", and "Split", corresponding to the risk of the event
	calculated by the model, and patient assignment to low versus
	high-risk groups, respectively.}
      \item{\code{dsname}:}{A character vector of length one; the name
	of the data set being modeled from a \code{MultiOmics} object.}
      \item{\code{SF}:}{Object of type \code{survfit}. See documentation
	for \code{\link{survfit}}.}
      \item{\code{riskModel}:}{Object of type \code{coxph}. See
	documentation for \code{\link{coxph}}.}
      \item{\code{splitModel}:}{Object of type \code{coxph}. See
	documentation for \code{\link{coxph}}.}

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
\value{

}
\author{Kevin R. Coombes \email{krc@silicovore.com}}
\note{
}
\seealso{
  \code{\link{getSizes}}
}
\examples{
data("TCGA-ESCA")
MO <- prepareMultiOmics(assemble, Outcome)
MO <- MO[c("ClinicalBin", "ClinicalCont", "RPPA"),]
set.seed(98765)
splitVec <- rbinom(nrow(Outcome), 1, 0.6)
trainD <- MO[, splitVec == 1]
testD <- MO[, splitVec == 0]

firstPass <- fitSingleModel(trainD, "RPPA",
                           "Days", "vital_status", "dead")
summary(firstPass)
plot(firstPass)
pre1 <- predict(firstPass, testD)
}
\keyword{classes}
