\name{matrix.utility}
\alias{flipud}
\alias{fliplr}

\title{OOMPA Matrix Utility Functions}
\description{
  Utility functions for manipulating matrices.
}
\usage{
flipud(x)
fliplr(x)
}
\arguments{
  \item{x}{a matrix}
}
\value{
  The \code{flipud} function returns a matrix the same size as \code{x},
  with the order of the rows reversed, so the matrix has been flipped
  vertically. The \code{fliplr} function returns a matrix the same size
  as \code{x} but flipped horizontally, with the order of the columns
  reversed. 
}
\author{Kevin R. Coombes <kcoombes@mdanderson.org>}
\examples{
mat <- matrix(1:6, 2, 3)
mat
flipud(mat)
fliplr(mat)
}
\keyword{array}
