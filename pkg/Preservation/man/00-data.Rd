\name{CAVA DAta}
\alias{prot5}
\alias{colvec}
\alias{paths-idex}
\alias{idex}
\alias{drs}
\docType{data}
\title{Dimension Reduction Applied to CAVA Daat Subset}
\description{
  A subset of protein expression data from a CITE-Seq experiment
  performed in the laboratory of Klaus Ley.
}
\usage{
data("prot5")
data("path-idex")
data("drs")
}
\format{
  The \code{prot5} data set contains two objects. one,also called
  \code{prot5} is a numerical matrix protein expression data having 51
  columns (proteins) and 3096 rows (individual cells). The second,
  \code{colvec} is a character vector holding coor assignments dividing
  the \code{prot5} data into four subclusters.

  The \code{paths-idex} set contains an integer vector, \code{idex} of
  length 301, indirectly defining a path by pointing, in order, to rows
  of the data marix.

  The \code{drs} set contains one object, \code{drs}, which is a list of
  four matrices, each with 2 columns and 3096 rows, representing the
  image of the 51-dimensioanl manifold after applying differnt dimension
  reduction algorithms,
}
\details{
  The \code{prot5} object is a subset of the protein data from the CAVA
  data set [Saigusa et al.]. The subset was slected because it represent
  a single coherent cluster of the CD4 T-cells in the experiment. We
  used hierarchical clustering with Ward's linkage rule and Euclidean
  distance to subdivide this cluster into four parts. The
  \code{colvec} vector records these cluster assignments. We then
  constructed a path in high dimendisonal space by performing principal
  components analysis and slecting all points close to the first
  principal component axis, ordered by their proection onto the
  axis. The (ordered) indices of the rows of these cells is sotred in
  the \code{idex} object. We also applid four dimension reduction
  algorithms (MDS and TSNE, using the \code{Rdimtools} package; UMAP,
  using the \code{umap} package, and DM, using the \code{diffusionMap}
  package) to the full \code{prot5} data set. The coordinates of the
  various rductions are stored inthe \code{drs} list.

}
\source{
  
}
\references{
  Saigusa et al., Sex Differences in Coronary Artery Disease and
  Diabetes Revealed by scRNA-Seq and CITE-Seq of Human CD4+ T
  Cells. Int J Mol Sci. 2022; 23(17):9875. doi: 10.3390/ijms23179875.
}
\seealso{
}
\examples{
data(prot5)
dim(prot5)
data(drs)
plot(drs[["MDS"]], col =colvec, pch = 16)
}
\keyword{datasets}
