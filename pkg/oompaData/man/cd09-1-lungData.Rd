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
  used to be available by FTP from the NIH, but can now only be found
  in the Gene Expression Omnibus at
  \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68571}.
  The data were log transformed by mapping the expression value x to
  log2(1+x). A subset of genes and of clinical annotation columns were
  selected to form this data set.
}
\references{
  Abrams ZB, Zucker M, Wang M, Asiaee Taheri A, Abruzzo LV, Coombes KR.\cr
  \emph{Thirty biologically interpretable clusters of transcription
    factors distinguish cancer type.}\cr
  BMC Genomics. 2018 Oct 11;19(1):738. doi: 10.1186/s12864-018-5093-z.

  Asiaee A, Abrams ZB, Nakayiza S, Sampath D, Coombes KR.\cr
  \emph{Explaining Gene Expression Using Twenty-One MicroRNAs.}\cr
  J Comput Biol. 2020 Jul;27(7):1157-1170. doi: 10.1089/cmb.2019.0321.
}
\keyword{datasets}

