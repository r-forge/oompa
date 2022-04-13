setClass("MultiplePLSCoxModels",
         slots = c(models = "list",
                   timevar = "character",
                   eventvar = "character",
                   eventvalue = "character"))

## As in MOFA, each data set must contain the same set of samples
validMultipleCoxModels <- function(object) {
  types <- all(sapply(object@models, function(MO) (
    inherits(MO, "SingleModel"))
  ))
  types
}
setValidity("MultiplePLSCoxModels", validMultipleCoxModels)

## summary method for MultipleCoxModels objects
setMethod("summary", "MultiplePLSCoxModels", function(object, ...) {
  cat("An object containing MultiplePLSCoxModels based on:\n",
      file = stdout())
  print(names(object@models))
})



## Need to define a class for the return value of this function
fitCoxModels <- function(object, timevar, eventvar, eventvalue, verbose = TRUE) {
  firstPass <- lapply(names(object@data), function(N) {
    if(verbose) cat(N, "\n", file = stderr())
    fitSingleModel(object, N, timevar, eventvar, eventvalue)
  })
  names(firstPass) <- names(object@data)
  new("MultiplePLSCoxModels",
      models = firstPass,
      timevar = timevar,
      eventvar = eventvar,
      eventvalue = eventvalue)
}

extendCoxModels <- function(object, firstPass, verbose = TRUE) {
  Components <- lapply(firstPass@models, function(DS) {
    DS@plsmod$tt
  })
  tempComps <- t(sapply(Components, dim))
  dimnames(tempComps) <- list(names(Components), c("nPats", "nComps"))
  sizing <- apply(tempComps, 2, sum)
  stores <- matrix(NA, nrow(object@outcome), sizing[2])
  rownames(stores) <- rownames(object@outcome)
  foo <- unlist(sapply(tempComps[,2], function(K) 1:K))
  colnames(stores) <- names(foo)
  ##
  for (N in names(Components)) {
    cat(N, "\n", file = stderr())
    X <- Components[[N]]
    colnames(X) <- paste(N, 1:ncol(X), sep = "")
    stores[rownames(X), colnames(X)] <- X
  }
  ##
  featgroups <- substring(colnames(stores), 1, -1 + nchar(colnames(stores)))
  wickedBig <- lapply(names(object@data), function(N) {
    if(verbose) cat(N, "\n", file = stderr())
    X <- object@data[[N]]
    allNA <- apply(X, 2, function(xcol) all(is.na(xcol)))
    X <- X[, !allNA]
    ident <- apply(X, 1, function(x) length(unique(x)))
    X <- t(X[ident > 1, ])
    Xout <- object@outcome[rownames(X),]
    Xstores <- stores[rownames(X),]
    mango <- lapply(names(object@data), function(M) {
      cat("\t", M, "\n", file = stderr())
      Y = Xstores[, featgroups == M]
      learn <- try(plsr(Y ~ X, 2))
      if (inherits(learn, "try-error")) return(NULL)
      extend <- predict(learn, X)
      list(learn = learn, extend = extend)
    })
    names(mango) <- names(object@data)
    t(sapply(mango, function(x) dim(x$extend)))
    slurp <- sapply(mango, function(x) x$extend[,,2])
    allPred <- do.call(cbind, slurp)
    list(mango = mango, allPred = allPred)
  })
  ##
  bilge <- sapply(wickedBig, function(x) x$allPred)
  class(bilge)
  sapply(bilge, dim)

  myArray <- array(NA, dim = c(nrow(stores), ncol(stores), length(bilge)))
  dimnames(myArray) <- list(rownames(stores), colnames(stores),
                            names(object@data))
  for (I in 1:length(bilge)) {
    B <- bilge[[I]]
    myArray[rownames(B), colnames(B), I] <- B
  }

  meanPreds <- apply(myArray, 1:2, mean, na.rm = TRUE)
}
