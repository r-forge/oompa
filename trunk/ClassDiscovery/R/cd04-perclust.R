################################################
setClass('PerturbationClusterTest',
         representation('ClusterTest',
                        f='function',
                        nTimes='numeric',
                        noise='numeric'))

PerturbationClusterTest <- function(data, FUN, nTimes=100, noise=1, verbose=TRUE, ...) {
  call <- match.call()
  if(inherits(data, 'ExpressionSet')) {
    data <- exprs(data)
  }
  N <- ncol(data)
  bootMatch <- matrix(0, nrow = N, ncol = N)
  for(i1 in 1:nTimes){
    tempData <- data + matrix(rnorm(N*nrow(data), 0, noise), ncol=N)
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
  new('PerturbationClusterTest', testResult,
      f=FUN, nTimes=nTimes, noise=noise)
}

setMethod('summary', signature(object='PerturbationClusterTest'),
          function(object, ...) {
  cat(paste('Number of perturbation samples: ', object@nTimes, '.\n', sep=''))
  cat(paste('Noise level: ', object@noise, '.\n', sep=''))
  callNextMethod()
})

