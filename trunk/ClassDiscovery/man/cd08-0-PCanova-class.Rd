\name{PCanova}
\alias{PCanova}
\alias{PCanova-class}
\alias{plot,PCanova,missing-method}
\alias{pltree,PCanova-method}
\alias{summary,PCanova-method}
\docType{class}
\title{Class "PCanova"}
\description{
  Implements the PCANOVA method for determining whether a putative group
  structure is truly reflected in multivariate data set.
}
\usage{
PCanova(data, classes, labels, colors, usecor=TRUE)

\S4method{plot}{PCanova,missing}(x, tag='', mscale=1, cex=1, \dots)
}
\arguments{
  \item{data}{either data frame or matrix with numeric values, or an
    \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} as defined
    in the BioConductor tools for analyzing microarray data.}
  \item{labels}{character vector used to label points in plots. The
    length of the \code{labels} vector should equal the number of
    columns (samples) in the \code{data} matrix. Since only the first
    character of each label is used in the plots, these should be unique.}
  \item{classes}{ A subset of the \code{labels} used to indicate
    distinct classes. Again, the method truncates each class indicator
    to a single letter.}
  \item{colors}{character vector containing color names; this should
    be the same length as the vector of \code{labels}.}
  \item{usecor}{logical scalar. If \code{TRUE}, standardize the rows of the
    data matrix before use.}
  \item{x}{object of class \code{PCanova}}
  \item{tag}{character string to name the object, used as part of the
    plot title.}
  \item{mscale}{A real number. This is a hack; for some reason, the
    projection of the sample vectors into the principal component space
    computed from the matrix of group means seems to be off by a factor
    approximately equal to the square root of the average number of
    samples per group. Until we sort out the correct formula, this term
    can be adjusted until the group means appear to be in the correct
    place in the plots.}
  \item{cex}{Character expansion factor used only in the plot legend on
    the plot of PC correlations.}
  \item{\dots}{additional graphical parameters passed on to \code{plot}
    when displaying the principal components plots.}
}
\details{
  The PCANOVA method was developed as part of the submission that won
  the award for best presentation at the 2001 conference on the Critical
  Assessment of Microarray Data Analysis (CAMDA;
  \url{http://www.camda.duke.edu}). The idea is to perform the
  equivalent of an analysis of variance (ANOVA) in principal component
  (PC) space. Let X(i,j) denote the jth column vector belonging to the
  ith group of samples. We can model this as X(i,j) = mu + tau(i) +
  E(i,j), where mu is the overall mean vector, tau(i) is the
  \dQuote{effects} vector for the ith group, and E(i,j) is the vector of
  residual errors. We can perform principal components analysis on the
  full matrix X containing all the columns X(i,j), on the matrix
  containing all the group mean vectors mu + tau(i), and on the residual
  matrix containing all the E(i,j) vectors. PCANOVA develops a measure
  (\dQuote{PC correlation}) for comparing these three sets of principal
  components. If the PC correlation is close to 1, then two principal
  component bases are close together; if the PC correlation is close to
  zero, then two principal components bases are dissimilar. Strong group
  structures are recognizable because the PC correlation between the
  total-matrix PC space and the group-means PC space is much larger than
  the PC correlation between the total-matrix PC space and the residual
  PC space. Weak or nonexistent group structures are recognizable
  because the relative sizes of the PC correlations is reversed.
}
\value{
  The \code{PCanova} function returns an object of the \code{PCanova} class.
}
\section{Objects from the Class}{
  Objects should be created by calling the \code{PCanova} generator function.
}
\section{Slots}{
  \describe{
    \item{\code{orig.pca}:}{A \code{matrix} containing the \code{scores}
    component from PCA performed on the total matrix.  All principal
    components analyses are performed using the \code{SamplePCA} class.}
    \item{\code{class.pca}:}{A \code{matrix} containing the
      \code{scores} component from PCA performed on the matrix of
      group-mean vectors. }
    \item{\code{resid.pca}:}{A \code{matrix} containing the
      \code{scores} component from PCA performed on the matrix of
      residuals. } 
    \item{\code{mixed.pca}:}{A \code{matrix} containing the projections
      of all the original vectoprs into the principal component space
      computed from the matrix of group mean vectors. }
    \item{\code{xc}:}{An object produced by performing hierarchical
      clustering on the total data matrix, using \code{hclust} with
      pearson distance and average linkage.  }
    \item{\code{hc}:}{An object produced by performing hierarchical
      clustering on the matrix of group means, using \code{hclust} with
      pearson distance and average linkage.}
    \item{\code{rc}:}{An object produced by performing hierarchical
      clustering on the matrix of residuals, using \code{hclust} with
      pearson distance and average linkage.}
    \item{\code{n}:}{An integer; the number of samples.}
    \item{\code{class2orig}:}{The \code{numeric} vector of PC
      correlations relating the total-matrix PCA to the group-means PCA.} 
    \item{\code{class2resid}:}{The \code{numeric} vector of PC
      correlations relating the residual PCA to the group-means PCA.}
    \item{\code{orig2resid}:}{The \code{numeric} vector of PC
      correlations relating the total-matrix PCA to the residual PCA.}
    \item{\code{labels}:}{A \code{character} vector of plot labels to
      indicate the group membership of samples.}
    \item{\code{classes}:}{A \code{character} vector of labels
      identifying the distinct groups. }
    \item{\code{colors}:}{A character vector of color names used to
      indicate the group membership fo samples in plots.}
    \item{\code{call}:}{An object of class \code{call} that records how
      the object was constructed.}
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = PCanova, y = missing)}: Plot the
      results of the PCANOVA test on the data. This uses \code{par} to
      set up a 2x2 layout of plots. The first three plots show the
      sample vectors (color-coded and labeled) in the space spanned by
      the first two principal components for each of the there PCAs. The
      final plot shows the three sets of PC correlations.  Colors in the
      first three plots are determined by the \code{colors} slot of the
      object, which was set when the object was created. Colors in the
      PC correlation plot are determined by the current values of
      \code{\link[oompaBase:oompaColor]{oompaColor$OBSERVED}},
      \code{\link[oompaBase:oompaColor]{oompaColor$EXPECTED}}, and
      \code{\link[oompaBase:oompaColor]{oompaColor$PERMTEST}}
    }
    \item{pltree}{\code{signature(x = PCanova)}: Produce dendrograms of
      the three hierarchical clusters of the samples, based on all the
      data, the group means, and the residuals. Since this method uses
      \code{par} to put all three dendrograms in the same window, it
      cannot be combined with other plots. }
    \item{summary}{\code{signature(object = PCanova)}: Write out a
      summary of the object. }
  }
}
\section{BUGS}{
  [1] The projection of the sample vectors into the principal component
  space of the group-means is off by a scale factor.  The \code{mscale}
  parameter provides a work-around.

  [2] The pltree method fails if you only supply two groups; this may be
  a failure in \code{hclust} if you only provide two objects to cluster.
}
\references{
  Examples of the output of PCANOVA applied to the NCI60 data set can be found
  at \url{http://bioinformatics.mdanderson.org/camda01.html}. The full
  description has not been published (out of laziness on the part of the
  author of this code). The only description that has appeared in print
  is an extremely brief description that can be found in the proceedings of
  the CAMDA 2001 conference.
}
\author{
  Kevin R. Coombes \email{krc@silicovore.com}
}
\seealso{
  \code{\link{SamplePCA}}
}
\examples{
showClass("PCanova")

## simulate data from three groups
d1 <- matrix(rnorm(100*10, rnorm(100, 0.5)), nrow=100, ncol=10, byrow=FALSE)
d2 <- matrix(rnorm(100*10, rnorm(100, 0.5)), nrow=100, ncol=10, byrow=FALSE)
d3 <- matrix(rnorm(100*10, rnorm(100, 0.5)), nrow=100, ncol=10, byrow=FALSE)
dd <- cbind(d1, d2, d3)
## colors that match the groups
cols <- rep(c('red', 'green', 'blue'), each=10)

## compute the PCanova object
pan <- PCanova(dd, c('red', 'green', 'blue'), cols, cols)
summary(pan)

## view the PC plots
plot(pan)

## view the dendrograms
pltree(pan, line=-0.5)

## compare the results when there is no underlying group structure
dd <- matrix(rnorm(100*50, rnorm(100, 0.5)), nrow=100, ncol=50, byrow=FALSE)
cols <- rep(c('red', 'green', 'blue', 'orange', 'cyan'), each=10)
pan <- PCanova(dd, unique(cols), cols, cols)
plot(pan, mscale=1/sqrt(10))

pltree(pan, line=-0.5)

## cleanup
rm(d1, d2, d3, dd, cols, pan)
}
\keyword{classes}
\keyword{multivariate}

