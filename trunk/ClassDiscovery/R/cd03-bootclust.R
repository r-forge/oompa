################################################
cutHclust <- function(data, k, method='average', metric = 'pearson') {
  tempCorr <- distanceMatrix(data, metric=metric)
  tempCorrTree <- hclust(tempCorr, method = method);
  tempCut <- cutree(tempCorrTree, k = k);
}

require(cluster)
cutPam <- function(data, k) {
  pamres <- pam(t(data), k)
  pamres$clustering
}

cutKmeans <- function(data, k) {
  kcent <- sample(ncol(data), k)
  kres <- kmeans(t(data), t(data[,kcent]))
  kres$cluster
}

repeatedKmeans <- function(data, k, nTimes) {
  n.samples <- ncol(data)
  kcent <- sample(n.samples, k)
  kres <- kmeans(t(data), t(data[,kcent]))
  withinss <- sum(kres$withinss)
  for (i in 1:nTimes) {
    tcent <- sample(n.samples, k)
    tres <- kmeans(t(data), t(data[,tcent]))
    if (sum(tres$withinss) < withinss) {
      kcent <- tcent
      kres <- tres
      withinss <- sum(kres$withinss)
    }
  }
  list(kmeans=kres, centers=kcent, withinss=withinss)
}

cutRepeatedKmeans <- function(data, k, nTimes) {
  x <- repeatedKmeans(data, k, nTimes)
  x$kmeans$cluster
}

################################################
setClass('BootstrapClusterTest',
         representation('ClusterTest',
                        f='function',
                        subsetSize='numeric',
                        nTimes='numeric'))


BootstrapClusterTest <- function(data, FUN, subsetSize, nTimes=100, verbose=TRUE, ...) {
  call <- match.call()
  if(inherits(data, 'ExpressionSet')) {
    data <- exprs(data)
  }
  N <- ncol(data)
  if (missing(subsetSize)) {
    subsetSize <- nrow(data)
  }
  subsetSize <- as.integer(subsetSize)
  bootMatch <- matrix(0, nrow = N, ncol = N)
  for(i1 in 1:nTimes){
    tempData <- data[sample(nrow(data), subsetSize, replace=TRUE),];
    if(verbose) {
      cat(paste('[', i1, '] ', nrow(tempData), ' ', sep=''))
      if (i1 %% 10 == 0) cat('\n')
    }
    tempCut <- FUN(tempData, ...)
    K <- max(tempCut)
    tempMatch <- matrix(0, nrow = N, ncol = N);
    for(i2 in 1:K){
      tempMatch[tempCut == i2, tempCut == i2] <- 1;
    }
    bootMatch <- bootMatch + tempMatch; 
  }
  dimnames(bootMatch) <- list(colnames(data), colnames(data))
  if(verbose) cat('\n')
  testResult <- new('ClusterTest', call = call, result = bootMatch/nTimes)
  new('BootstrapClusterTest', testResult,
      f=FUN, nTimes=nTimes, subsetSize=subsetSize)
}

setMethod('summary', 'BootstrapClusterTest', function(object, ...) {
  cat(paste('Number of bootstrap samples: ', object@nTimes, '.\n', sep=''))
  cat(paste('Number of rows sampled: ', object@subsetSize, '.\n', sep=''))
  callNextMethod()
})





