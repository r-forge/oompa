\name{BumSummary-class}
\docType{class}
\alias{BumSummary}
\alias{BumSummary-class}
\alias{show,BumSummary-method}

\title{The BumSummary class }
\description{
  An implementation class. Users are not expected to create these objects directly;
  they are produced as return objects from the summary method for \code{Bum}.
}
\section{Slots}{
  \describe{
    \item{\code{bum}:}{Object of class \code{Bum} }
    \item{\code{estimates}:}{Object of class \code{data.frame}}
    \item{\code{Fhat}:}{Object of class \code{numeric}}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "BumSummary")}: Print the
      object, which contains a summary of the underlying \code{Bum}
      object. The summary contains a data frame with estimates of the
      fraction of true positives (TP), false postives (FP), true negatives
      (TN) and false negatives (FN) at the set of p-value cutoffs
      specified in the call to the summary method.}
  }
}
\author{Kevin R. Coombes <kcoombes@mdanderson.org>}
\seealso{
 \code{\link{Bum}}
}
\examples{
# See the examples for the Bum class.
}
\keyword{classes}
