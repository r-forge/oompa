\name{fitCoxModel-class}
\docType{class}
\alias{fitCoxModel}
\alias{fitCoxModel-class}
\alias{predict,fitCoxModel-method}
\alias{summary,fitCoxModel-method}
\alias{plot,fitCoxModel-method}
\alias{fitSingleModel}
\alias{getSizes}

\title{Class \code{"fitCoxModel"}}
\description{
  The \code{fitCoxModel} object class ...
  The \code{fitCoxModels} function takes in an object of \code{MultiOmics} class and returns a new object of \code{MultiplePLSCoxModels} class.
}
\section{Objects from the Class}{

}
\usage{
fitCoxModels(object, timevar, eventvar, eventvalue, verbose = TRUE)
\S4method{summary}{fitCoxModel}(object, \dots)
\S4method{plot}{fitCoxModel,missing}(x, y, \dots)
}
\arguments{
  \item{object}{an object of class \code{MultiOmics} for fitting the model.}
  \item{timevar}{a column in the \code{MultiOmics} object in the \code{outcome} dataframe containing the time-to-event.}
  \item{eventvar}{a column in the \code{MultiOmics} object in the \code{outcome} dataframe containing the event.}
  \item{eventvalue}{a character string specifying the value of the event.}
  \item{verbose}{prints messages describing the model fitting process. Default is TRUE.}
  
}
\section{Slots}{
  \describe{
      \item{\code{plsmod}:}{Object of class \code{plsRcoxmodel} containing the fitted model.}
      \item{\code{Xout}:}{Object of type \code{data.frame} containing the original \code{outcome} dataframe and additional columns for "Risk", and "Split", corresponding to the risk of the event calculated by the model, and patient assignment to low versus high-risk groups, respectively.}
      \item{\code{SF}:}{Object of type \code{survfit}. See documentation for \code{\link{survfit}}.}
      \item{\code{riskModel}:}{Object of type \code{coxph}. See documentation for \code{\link{coxph}}.}
      \item{\code{splitModel}:}{Object of type \code{coxph}. See documentation for \code{\link{coxph}}.}

  }
}
\section{Methods}{
  \describe{
      \item{\code{\link{plot}}:}{Plots Kaplan-Meier curves for each omics dataset split into Low Risk and High Risk groups.}
      \item{\code{\link{summary}}:}{Returns a description of the \code{MultiplePLSCoxModels} object and the names of the omics datasets used to build the model.}
  
      \item{\code{\link{predict}}:}{not sure what to put here, wasn't able to run predict in the vignette.}
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

set.seed(98765)
splitVec <- rbinom(nrow(Outcome), 1, 0.6)
trainD <- MO[, splitVec == 1]
testD <- MO[, splitVec == 0]

firstPass <- fitCoxModels(trainD, timevar = "Days",
                          eventvar = "vital_status", eventvalue = "dead")
summary(firstPass)
plot(firstPass)
#pre1 <- predict(firstPass, testD)
}
\keyword{classes}
