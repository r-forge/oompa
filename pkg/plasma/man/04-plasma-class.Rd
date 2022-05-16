\name{plasma-class}
\docType{class}
\alias{plasma-class}
\alias{plot,plasma,missing-method}
\alias{plasma}

\title{Class \code{"plasma"}}
\description{
  The \code{plasma} object class is returned after running the \code{plasma} function.
  The \code{plasma} function uses the \code{PLSRCox} components from one dataset as the predictor variables and the features of another dataset as the response variables to fit a partial least squares regression (\code{plsr}) model. Then, we take the average of the predictions to create a final matrix of samples versus components.

  In practice, the user will feed in a MultiOmics object and a MultiplePLSCoxModels object to the \code{plasma} function. The matrix of components described earlier is then used to fit a Cox Proportional Hazards (\code{coxph}) model with AIC stepwise variable selection to return a final object of class \code{plasma} which includes a \code{coxph} model with a reduced number of predictors.
}
\section{Objects from the Class}{
Objects should be defined using the \code{plasma} function.
}
\usage{
plasma(object, multi)
\S4method{plot}{plasma}(x, y, \dots)
\S4method{predict}{plasma}(x, y, \dots)
\S4method{heatmap}{plasma}(x, y, \dots)
}
\arguments{
  \item{object}{an object of the \code{MultiOmics} class.}
  \item{multi}{an object of the \code{MultiplePLSCoxModels} class.}  
  \item{meanPredictions}{an matrix of meanPredictions stored in an object of the \code{plasma} class.}
  \item{\dots}{parameters for \code{plot} and \code{predict} are inherited from the \code{plasmaPredictions} class.}
}
\section{Slots}{
  \describe{
      \item{\code{traindata}:}{An object of class \code{MultiOmics} used for training the model.}
      \item{\code{compModels}:}{A list containing objects in the form of \code{plsr}.}
      \item{\code{fullModel}:}{A coxph object with variables (components) selected via AIC stepwise selection.}
  }
}
\section{Methods}{
  \describe{
      \item{\code{\link{plot}}:}{Plots a Kaplan-Meier curve of the final \code{coxph} model that has been categorized into \dQuote{low risk} and \dQuote{high risk} based whether it is higher or lower, respectively, than the median value of risk.}
      \item{\code{\link{predict}}:}{spaceholder text}
      \item{\code{\link{heatmap}}:}{spaceholder text}
  }
}
\value{
}
\author{Kevin R. Coombes \email{krc@silicovore.com}}
\note{
}
\seealso{plasmaPredictions, plsr}
\examples{
}
\keyword{classes}
