# Copyright (C) Kevin R. Coombes, 2007-2012

# pcanova

#########################################################
# Create yet another class of objects

setClass('PCanova',
         representation(orig.pca  = 'matrix',
                        class.pca = 'matrix',
                        mixed.pca = 'matrix',
                        resid.pca = 'matrix',
                        xc = 'hclust',
                        hc = 'hclust',
                        rc = 'hclust',
                        n = 'numeric',
                        class2orig  = 'numeric',
                        class2resid = 'numeric',
                        orig2resid  = 'numeric',
                        labels  = 'character',
                        classes = 'character',
                        colors  = 'character',
                        call='call'))

setMethod('summary', signature(object='PCanova'),
          function(object, ...) {
  cat(paste('An object of the', class(object),  'class with', object@n, 'groups.\n'))
  cat(paste('Call:\n\t', as.character(list(object@call)), '\n', sep=''))
})

PCanova <- function(data, classes, labels, colors, usecor=TRUE) {
  # note 1: classes must be a subset of labels

  call <- match.call()
  if(inherits(data, 'ExpressionSet')) {
    data <- exprs(data)
  }
  data <- as.matrix(as.data.frame(data))
  labels <- substring(labels, 1, 1)
  classes <- unique(substring(classes, 1, 1))
  # function that yields column vectors with all ones
  ones <- function(rows, cols=1) { matrix(1, rows, cols) }

  n <- dim(data)[2]	# number of samples
  g <- length(classes)	# number of groups
  p <- dim(data)[1]	# number of variables
  qq <- matrix(0, g, n)	# transpose of n x g class indicator matrix
  for(i in 1:n) {
    qq[classes==labels[i], i] <- 1
  }
  nn <- qq %*% ones(n) 			# number of samples in each class
  xsums <- data %*% t(qq)			# matrix of class sums
  cc <- xsums %*% diag(as.vector(1/nn))	# matrix of class means
  xm <- data %*% ones(n)/n			# vector of mean values
  zz <- data - (xm %*% t(ones(n))) 	# center the rows
  vv <- (cc - (xm %*% t(ones(g))))	# matrix of centered class means
  ss <- cc %*% qq - xm %*% t(ones(n))	#
  ss0 <- cc - xm %*% t(ones(g))
  rr <- data  - (cc %*% qq)		# matrix of residuals

  # first step: do PCA on original data
  xd <- SamplePCA(zz, usecor=usecor, center=FALSE)
  orig.pca <- xd@scores

  # second step: form the vectors of class medians into a matrix
  yd <- SamplePCA(ss0, usecor=usecor, center=FALSE)
  class.pca <- yd@scores
  mixed.pca <- t(zz) %*% yd@components

  # third step: compute the residual matrix
  rd <- SamplePCA(rr, usecor=usecor, center=FALSE)
  resid.pca <- rd@scores

  # fourth step: compute transformation matrices from one PC to another
  # note: we only keep a number of PCS equal to the number of classes
  n <- length(classes)
  trans.class.orig <-  t(t(yd@components) %*% xd@components)[1:n, 1:n]
  trans.class.resid <- t(t(yd@components) %*% rd@components)[1:n, 1:n]
  trans.orig.resid <-  t(t(xd@components) %*% rd@components)[1:n, 1:n]

  class2orig <- rep(0, n)
  class2resid <- rep(0, n)
  orig2resid <- rep(0, n)
  for (i in 1:n) {
    class2orig[i]  <- 1-acos((sqrt(sum(trans.class.orig[i, 1:i]^2)) +
                              sqrt(sum(trans.class.orig[1:i,i]^2)))/2)  /(pi/2)
    class2resid[i] <- 1-acos((sqrt(sum(trans.class.resid[i, 1:i]^2)) +
                              sqrt(sum(trans.class.resid[1:i,i]^2)))/2)/(pi/2)
    orig2resid[i]  <- 1-acos((sqrt(sum(trans.orig.resid[i, 1:i]^2)) +
                              sqrt(sum(trans.orig.resid[1:i,i]^2)))/2)  /(pi/2)
  }
  # finally do some clustering on correlation
  hc <- hclust(distanceMatrix(cc, 'pearson'), method='average')
  xc <- hclust(distanceMatrix(data, 'pearson'), method='average')
  rc <- hclust(distanceMatrix(rr, 'pearson'), method='average')
  new('PCanova', orig.pca = orig.pca, class.pca = class.pca, mixed.pca = mixed.pca,
      resid.pca = resid.pca, xc = xc, hc = hc, rc = rc, n = n,
      class2orig = class2orig, class2resid = class2resid, orig2resid = orig2resid,
      labels = labels, classes = classes, colors = colors, call=call)
}

setMethod('plot', signature('PCanova', 'missing'),
          function(x, tag='', mscale=1, cex=1, ...) {
  opar <- par(mfrow=c(2,2))
  on.exit(par(opar))
  plot(x@orig.pca[,1], x@orig.pca[,2], pch=x@labels, col=x@colors,
                xlab='comp. 1', ylab='comp. 2', ...)
  title('PCA by Element')
  groups <- sort(unique(x@labels))
  reps <- unlist(lapply(groups, function(g, lab) {sort(grep(g, lab))[1]}, x@classes))
  plot(c(x@mixed.pca[,1]/mscale, x@class.pca[,1]),
       c(x@mixed.pca[,2]/mscale, x@class.pca[,2]), type='n',
       xlab='comp. 1', ylab='comp. 2', ...)
  points(x@mixed.pca[,1]/mscale, x@mixed.pca[,2]/mscale,
         pch=x@labels, col=x@colors)
  points(x@class.pca[reps,1], x@class.pca[reps,2],
         pch=groups, col=rep('black', length(groups)), cex=2)
  title('PCA by Class')
  plot(x@resid.pca[,1], x@resid.pca[,2], pch=x@labels, col=x@colors,
       xlab='comp. 1', ylab='comp. 2', ...)
  title('Residual PCA')
  plot(c(0, x@n), c(0,2), type='n', ylim=c(0, 1),
       xlab='Number of Components', ylab='PC Correlation')
  points(cumsum(x@class2resid)/(1:x@n), type='b', col=oompaColor$OBSERVED, pch=15)
  points(cumsum(x@class2orig)/(1:x@n), type='b', col=oompaColor$PERMTEST, pch=16)
  points(cumsum(x@orig2resid)/(1:x@n), type='b', col=oompaColor$EXPECTED, pch=17)
  abline(h=1/2)
  legend((x@n-1)/2, 1, c('d(Class, Residual)', 'd(Total, Class)', 'd(Total, Residual)'),
         col=c(oompaColor$OBSERVED, oompaColor$PERMTEST, oompaColor$EXPECTED),
         pch=c(15, 16, 17), cex=cex)
  title(paste('Cumulative PC-ANOVA', tag))
  invisible(x)
})

setMethod('pltree', signature(x='PCanova'),
          function(x, ...) {
  opar <- par(mfrow=c(3,1))
  on.exit(par(opar))
  plotColoredClusters(x@xc, labs=x@labels, cols=x@colors, ...)
  title('Clustering the Samples')
  plotColoredClusters(x@rc, labs=x@labels, cols=x@colors, ...)
  title('Clustering the Residuals')
  plotColoredClusters(x@hc, labs=x@classes, cols=x@colors, ...)
  title('Clustering the Classes')
  invisible(x)
})

