\name{lungData-data}
\alias{lungData-data}
\alias{lungData}
\alias{lung.clinical}
\alias{lung.dataset}
\docType{data}
\title{Lung Cancer Gene Expression Dataset}
\description{
  This data set contains clinical annotations and the log expression of
  150 genes for a set of 444 lung cancer patients.  The 150 genes were
  selected randomly from a larger Affymetrix U133A dataset. 
}
\usage{
data(lungData)
}
\format{A data matrix (\code{lung.dataset}) containing the log
  expression of 150 genes (rows) in 444 lung tumor samples (columns),
  along with a data frame (lung.clinical) containing clinical
  annotations of the patients.} 
\source{
  Supporting data for the Nature Medicine paper by Shedden et al. was
  downloaded from the (now defunct) caArray web site. The original data
  are still available by FTP from 
  \url{ftp://caftpd.nci.nih.gov/pub/caARRAY/experiments/caArray_beer-00153/}
  and in the Gene Expression Omnibus at
  \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68571}.
  The data were log transformed by mapping the expression value x to
  log2(1+x). A subset of genes and of clinical annotation columns were
  selected to form this data set.
}
\references{
  Shedden K, Taylor JM, et al.\cr
  \emph{Gene expression-based survival prediction in lung adenocarcinoma:
  a multi-site, blinded validation study.}\cr
  (Director's Challenge Consortium for the Molecular Classification of
  Lung Adenocarcinoma).\cr
  Nature Medicine 2008; 14(8):822-7.
}
\keyword{datasets}

