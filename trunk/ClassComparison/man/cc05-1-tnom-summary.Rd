\name{TNoMSummary-class}
\docType{class}
\alias{TNoMSummary}
\alias{TNoMSummary-class}
\alias{show,TNoMSummary-method}

\title{The TNoMSummary Class}
\description{
  An implementation class. Users are not expected to create these objects directly;
  they are produced as return objects from the summary method for \code{TNoM}.
}
\section{Slots}{
  \describe{
    \item{\code{TNoM}:}{Object of class \code{TNoM} ~~ }
    \item{\code{counts}:}{Object of class \code{numeric} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = TNoMSummary)}:  Print the
      object, which contains a summary of the underlying \code{TNoM}
      object. In particular, the summary reports the number of genes
      acheiving each possible number of misclassifications.}
  }
}
\author{Kevin R. Coombes <kcoombes@mdanderson.org>}
\seealso{
 \code{\link{TNoM}}
}
\examples{
# See the examples for the TNoM class.
}
\keyword{classes}
