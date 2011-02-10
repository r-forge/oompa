# adds additional (common microarray) metrics to 'dist'
distanceMatrix <- function(dataset, metric, ...) {
  if(inherits(dataset, "ExpressionSet")) {
    dataset <- exprs(dataset)
  }
  METRICS <- c('pearson', 'sqrt pearson', 'spearman', 'weird',
               'absolute pearson', 'uncentered correlation')
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
    distance <- switch(metric,
                       pearson = as.dist((1-cor(dataset))/2),
                       'sqrt pearson' = as.dist(sqrt(1-cor(dataset))),
                       weird = dist(cor(dataset)),
                       'absolute pearson' = as.dist((1-abs(cor(dataset)))),
                       spearman = as.dist((1-cor(apply(dataset, 2, rank)))/2),
                       'uncentered correlation' = uncent(dataset))
  }
  distance
}

# Input is typically a data frame of log-transformed normalized
# expression data. Produces three hierarchical clusters, one from
# Euclidean distance, one from correlation, and one from manhattan
# distance of on-or-off genes.
cluster3 <- function(data, eps=logb(1,2), name='', labels=dimnames(data)[[2]]) {
  if(inherits(data, 'ExpressionSet')) {
    d <- exprs(data)
  } else {
    d <- data
  }
  plot(c(0,4), c(0,3), xaxt='n', yaxt='n', type='n', xlab='', ylab='')
  text(2, 2, name, cex=1.2)
  text(2, 1, 'Cluster Analysis', cex=1.2)
  plclust(hclust(distanceMatrix(d, 'euclid')), labels=labels)
  title(paste('Clustering on Euclidean distance'))
  plclust(hclust(distanceMatrix(d, 'pearson')), labels=labels)
  title('Clustering on correlation')
  plclust(hclust(distanceMatrix(d > eps, 'manhattan')), labels=labels)
  title('Clustering on presence or absence of genes')
  invisible(data)
}

plotColoredClusters <- function(hd, labs, cols, cex=0.8, main='', line=0, ...)
{
  plot(hd, hang = -1, labels = rep("", length(labs)), cex=cex, main=main, ...)
  mtext(labs, side=1, line=line, at=order(hd$order), col=cols, las=2, cex=cex)
}

