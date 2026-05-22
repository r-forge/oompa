\name{xkcd}
\alias{xkcd}
\docType{data}
\title{Color Names From the Color Survey by Randall Munroe of XKCD}
\description{
  A data frame mapping hex codes for 949 colors to their consensus
  names in the XKCD color survey.
}
\usage{data(xkcd)}
\format{
  A data frame with two columns and 949 rows.
}
\details{
  This data set contains the names and hex codes for 949 color names
  identified as consistent in the Color Survey run by Randall Munroe at
  XKCD.  During the online survey, more than five million color
  swatches were named across 222,500 user sessions. In the blog post
  describing the results, Munroe describes the method he used to come up
  with a list of "consensus" color names. He also says the text file
  should contain 954 color names, but the downloaded text file only
  contained 949. Data were reformatted to be usable as a data frame in
  R, with underscores inserted in place of spaces in multi-word color
  names.
}
\source{
  Our main source was the text file linked to from
  \url{https://blog.xkcd.com/2010/05/03/color-survey-results/}.
}
\seealso{
  \code{\link{xkcdNames}}
}
\examples{
data(xkcd)
head(xkcd)
}
\keyword{datasets}
\keyword{color}
