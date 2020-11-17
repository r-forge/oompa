# Copyright (C) Kevin R. Coombes, 2007-2012

unitize <- function(mat) {
  enorm <- sqrt(apply(mat^2, 2, sum))
  sweep(mat, 2, enorm, "/")
}

# adds additional (common microarray) metrics to 'dist'
distanceMatrix <- function(dataset, metric, ...) {
  if(inherits(dataset, "ExpressionSet")) {
    dataset <- Biobase::exprs(dataset)
  }
  METRICS <- c('pearson', 'sqrt pearson', 'spearman', 'weird',
               'absolute pearson', 'uncentered correlation',
               'cosine')
  m <- pmatch(metric, METRICS, nomatch=NA)
  if (is.na(m)) {
    distance <- dist(t(dataset), method=metric, ...)
  } else {
    metric <- METRICS[m]
    uncent <- function(dataset) {
      temp <- sqrt(apply(dataset^2, 2, mean))
      dataset <- sweep(dataset, 2, temp, '/')
      temp <- t(dataset) %*% dataset/nrow(dataset)
      as.dist(1 - temp)
    }
    cosine <- function(dataset) {
      uni <- unitize(dataset) # unit norm
      cosd <- t(uni) %*% uni        # similarity
      cosineDistance <- as.dist(1 - cosd)
    }
    distance <- switch(metric,
                       pearson = as.dist((1-cor(dataset))/2),
                       'sqrt pearson' = as.dist(sqrt(1-cor(dataset))),
                       weird = dist(cor(dataset)),
                       'absolute pearson' = as.dist((1-abs(cor(dataset)))),
                       spearman = as.dist((1-cor(apply(dataset, 2, rank)))/2),
                       'uncentered correlation' = uncent(dataset),
                       cosine = cosine(dataset))
  }
  distance
}

# Input is typically a data frame of log-transformed normalized
# expression data. Produces three hierarchical clusters, one from
# Euclidean distance, one from correlation, and one from manhattan
# distance of on-or-off genes.
cluster3 <- function(data, eps=logb(1,2), name='', labels=dimnames(data)[[2]]) {
  if(inherits(data, 'ExpressionSet')) {
    d <- Biobase::exprs(data)
  } else {
    d <- data
  }
  plot(c(0,4), c(0,3), xaxt='n', yaxt='n', type='n', xlab='', ylab='')
  text(2, 2, name, cex=1.2)
  text(2, 1, 'Cluster Analysis', cex=1.2)
  plot(hclust(distanceMatrix(d, 'euclid')), labels=labels, main="")
  title(paste('Clustering on Euclidean distance'))
  plot(hclust(distanceMatrix(d, 'pearson')), labels=labels, main="")
  title('Clustering on correlation')
  plot(hclust(distanceMatrix(d > eps, 'manhattan')), labels=labels, main="")
  title('Clustering on presence or absence of genes')
  invisible(data)
}

plotColoredClusters <- function(hd, labs, cols, cex=0.8, main='', line=0, ...)
{
  plot(hd, hang = -1, labels = rep("", length(labs)), cex=cex, main=main, ...)
  mtext(labs, side=1, line=line, at=order(hd$order), col=cols, las=2, cex=cex)
}



pcc <- function(hd, colors=NULL, ...) {
  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))
  if (!is.null(colors)) {
    L <- length(colors)
    lwid <- c(6, 1)
    lhei <- c(6, rep(0.5, L))
    lmat <- matrix(c(L+1, 1:L, L+2, rep(0,L)), ncol=2)
    layout(lmat, widths=lwid, heights=lhei, respect=FALSE)
    for (i in 1:L) {
      par(mar = c(0.5, 0.5, 0, 0))
      x <- colors[[i]]
      image(matrix(as.numeric(x$fac)[hd$order], ncol=1), col = x$col,
            axes = FALSE, xaxt='n', yaxt='n', xlab='', ylab='')
      mtext(side = 4, names(colors)[i], line = 0.5, las=1)
    }
    par(mar=c(5, 0.5, 7, 0))
  }
  plot(as.dendrogram(hd), yaxt='n', xaxs='i', ...)
  if (!is.null(colors)) {
    L <- length(colors)
    N <- sum(J <- unlist(lapply(colors, function(x) length(levels(x$fac)))))
    plot(c(0,1), c(1, N+L+1), type='n', xlab='', xaxt='n', ylab='', yaxt='n')
    for (i in 1:length(colors)) {
      y <- i + cumsum(J)[i]
      text(0.1, y+0.5, names(J)[i], adj=0)
      legend(0.3, y, levels(colors[[i]]$fac), col=colors[[i]]$col, pch=15)
    }
  } else {
    plot(c(0,1), c(0,1), type='n', xlab='', xaxt='n', ylab='', yaxt='n')
  }
  layout(matrix(1), widths=1, heights=1)
  invisible(hd)
}

if(FALSE) {
  fakedata <- matrix(rnorm(200*30), ncol=30)
  colnames(fakedata) <- paste("P", 1:30, sep='')

  faccol <- list(fac=factor(rep(c("A", "B"), each=15)),
                 col=c(A='red', B='green'))
  fac2col <- list(fac=factor(rep(c("X", "Y", "Z"), times=10)),
                 col=c(X='cyan', Y='orange', Z='purple'))
  BA <- faccol
  BA$col <- c(A='blue', B='yellow')
  colors <- list(AB=faccol, XYZ=fac2col, "tricky long name"=fac2col,
                 another=BA)

  hc <- hclust(distanceMatrix(fakedata, "pearson"), "ward")

  pcc(hc, colors)
}
