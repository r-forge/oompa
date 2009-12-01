\name{MultiLinearModel-class}
\docType{class}
\alias{MultiLinearModel}
\alias{multiTukey}
\alias{MultiLinearModel-class}
\alias{anova,MultiLinearModel-method}
\alias{hist,MultiLinearModel-method}
\alias{plot,MultiLinearModel,missing-method}
\alias{plot,MultiLinearModel,ANY-method}
\alias{summary,MultiLinearModel-method}

\title{The MultiLinearModel Class }
\description{
  Class to fit multiple (row-by-row) linear (fixed-effects) models on
  microarray or proteomics data.
}
\usage{
MultiLinearModel(form, clindata, arraydata)
\S4method{summary}{MultiLinearModel}(object, ...)
\S4method{hist}{MultiLinearModel}(x, xlab='F Statistics', main=NULL, ...)
\S4method{plot}{MultiLinearModel,missing}(x, y, ylab='F Statistics', ...)
\S4method{plot}{MultiLinearModel,ANY}(x, y, xlab='F Statistics',
 ylab=deparse(substitute(y)), ...)
\S4method{anova}{MultiLinearModel}(object, ob2, ...)
multiTukey(object, alpha)
}
\arguments{
  \item{form}{ A \code{formula} object specifying the linear model}
  \item{clindata}{ Either a data frame of "clinical" or other
    covariates or an \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}. }
  \item{arraydata}{ A matrix or data frame of values to be explained by
    the model.  If \code{clindata} is an \code{ExpressionSet}, then
    \code{arraydata} can be omitted, since it is assumed to be part of
    the \code{ExpressionSet}.}
  \item{object}{A \code{MultiLinearModel} object}
  \item{ob2}{Another \code{MultiLinearModel} object}
  \item{x}{A \code{MultiLinearModel} object}
  \item{y}{An optional numeric vector}
  \item{xlab}{Label for the x-axis}
  \item{ylab}{Label for the y-axis}
  \item{main}{Graph title}
  \item{\dots}{Optional graphical or other parameters to generic functions}
  \item{alpha}{A real number between 0 and 1; the significance level for
  the Tukey test.}
}
\section{Creating Objects}{
  Objects should be created by calling the \code{MultiLinearModel}
  function.  The first argument is a \code{\link{formula}} specifying
  the linear model, in the same manner that it would be passed to
  \code{\link{lm}}.  We will fit the linear model separately for each
  row in the \code{arraydata} matrix.  Rows of \code{arraydata} are
  attached to the \code{clindata} data frame and are always referred to
  as "Y" in the formulas. In particular, this implies that
  \code{clindata} can not include a column already called "Y".  Further,
  the implementation only works if "Y" is the response variable in the
  model.
}
\section{Multiple linear models with "ExpressionSet" objects}{
  The BioConductor packages uses an \code{ExpressionSet} to combine microarray
  data and clinical covariates (known in their context as
  \code{\link[Biobase]{phenoData}} objects) into a single structure.
  You can call \code{MultiLinearModel} using an \code{ExpressionSet} object
  for the \code{clindata} argument.  In this case, the function extracts
  the \code{phenoData} slot of the \code{ExpressionSet} to use for the
  clinical covariates, and extracts the \code{exprs} slot of the
  \code{ExpressionSet} object to use for the array data.
}
\section{Slots}{
  \describe{
    \item{\code{call}:}{A \code{call} object describing how the object
      was constructed. }
    \item{\code{model}:}{The \code{formula} object specifying the linear
    model. }
    \item{\code{F.statistics}:}{A numeric vector of F-statistics
      comparing the linear model to the null model. }
    \item{\code{p.values}:}{A numeric vector containing the p-values
      associated to the F-statistics.}
    \item{\code{coefficients}:}{A \code{matrix} of the coefficients in
      the linear models. }
    \item{\code{predictions}:}{A \code{matrix} of the (Y-hat) values
      predicted by the models. }
    \item{\code{sse}:}{A numeric vector of the sum of squared error
      terms from fitting the models. }
    \item{\code{ssr}:}{A numeric vector of the sum of squared
      regression terms from the model. }
    \item{\code{df}:}{A numeric vector (of length two) containing the
      degrees of freedom for the F-tests.}
  }
}
\section{Methods}{
  \describe{
    \item{summary(object, ...)}{Write out a summary of the object.}
    \item{hist(x, xlab='F Statistics', main=NULL, ...)}{Create a
      histogram of the F-statistics.}
    \item{plot(x, ylab='F Statistics', ...)}{Plot the F-statistics as a
      function of the row index.}
    \item{plot(x, y, xlab='F Statistics', ylab=deparse(substitute(y)),
      ...)}{Plot the F-statistics against the numeric vector \code{y}.}
    \item{anova(object, ob2, ...)}{Perform row-by-row F-tests comparing
      two different linear models.}
  }
}
\section{Details}{
  The \code{MultiLinearModel} constructor computes row-by-row F-tests
  comparing each linear model to the null model Y ~ 1.  In many
  instances, one wishes to use an F-test to compare two different linear
  models. For instance, many standard applications of analysis of
  variance (ANOVA) can be described using sucha compoarison between two
  different linear models. The \code{anova} method for the
  \code{MultiLinearModel} class performs row-by-row F-tests comparing
  two competing linear models.
  
  The implementation of \code{MultiLinearModel} does not take the naive
  approach of using either \code{\link{apply}} or a
  \code{\link{for}}-loop to attach rows one at a time and fit separate
  linear models. All the models are actually fit simultaneously by a
  series of matrix operations, which greatly reduces the amount of time
  needed to compute the models. The constraint on the column names in
  \code{clindata} still holds, since one row is attached to allow
  \code{\link{model.matrix}} to determine the contrasts matrix.
}
\value{
  The \code{anova} method returns a data frame.  The rows in the data
  frame correpsonds to the rows in the \code{arraydata} object that was
  used to construct the \code{MultiLinearModel} objects.  The first
  column contains the F-statistics and the second column contains the
  p-values.

  The \code{multiTukey} function returns a vector whosem length equals
  the number of rows in the \code{arraydata} object used to construct
  the \code{MultiLinearModel}. Assuming that the overall F-test was
  significant, differences in group means (in each data row) larger than
  this value are significant by Tukey's test for honestly significant
  difference. (Of course, that statement is incorrect, since we haven't
  fully corrected for multiple testing. Our standard practice is to take
  the p-values from the row-by-row F-tests and evaluate them using the
  beta-uniform mixture model (see \code{\link{Bum}}). For the rows that
  correspond to models whose p-values are smaller than the \code{Bum}
  cutoff, we simply use the Tukey HSD values without further
  modification.)
}
\author{Kevin R. Coombes <kcoombes@mdanderson.org>}
\seealso{
  \code{\link{MultiTtest}}, \code{\link{MultiWilcoxonTest}},
  \code{\link{Bum}}, \code{\link{lm}}, \code{\link{anova}}.
}
\examples{
ng <- 10000
ns <- 50
dat <- matrix(rnorm(ng*ns), ncol=ns)
cla <- factor(rep(c('A', 'B'), 25))
cla2 <- factor(rep(c('X', 'Y', 'Z'), times=c(15, 20, 15)))
covars <- data.frame(Grade=cla, Stage=cla2)
res <- MultiLinearModel(Y ~ Grade + Stage, covars, dat)
summary(res)
hist(res, breaks=101)
plot(res)
plot(res, res@p.values)

graded <- MultiLinearModel(Y ~ Grade, covars, dat)
summary(graded)

hist(graded@p.values, breaks=101)
hist(res@p.values, breaks=101)

oop <- anova(res, graded)
hist(oop$p.values, breaks=101)

# cleanup
rm(ng, ns, dat, cla, cla2, covars, res, graded, oop)
}
\keyword{classes}
\keyword{multivariate}
\keyword{htest}
