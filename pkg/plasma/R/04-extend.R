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

setClass("plasma",
         slots = c(meanPredictions = "matrix",
                   riskDF = "data.frame",
                   fullModel = "coxph",
                   riskModel = "coxph",
                   splitModel = "coxph",
                   SF = "survfit"))

## object = MultiOmics
## multi = MultiplePLSCoxModels
plasma <- function(object, multi) {
  meanPreds <- extendCoxModels(object, multi)
  colms <- c(multi@timevar, multi@eventvar)
  riskDF <- data.frame(object@outcome[, colms], meanPreds)
  form <- formula(paste("Surv(", multi@timevar, ",", multi@eventvar, "== \"",
                        multi@eventvalue, "\") ~ .", sep = ""))
  model <- coxph(form, data = riskDF)
  model1 <- step(model, trace = 0)
  riskDF$Risk <- predict(model1)
  lh <- c("Low", "High")
  riskDF$Split <- factor(lh[1 + 1*(riskDF$Risk >= median(riskDF$Risk))],
                          levels = lh)
  form <- formula(paste("Surv(", multi@timevar, ",", multi@eventvar, "== \"",
                        multi@eventvalue, "\") ~ Risk", sep = ""))
  rmodel <- coxph(form, riskDF)
  form <- formula(paste("Surv(", multi@timevar, ",", multi@eventvar, "== \"",
                        multi@eventvalue, "\") ~ Split", sep = ""))
  smodel <- coxph(form, riskDF)
  SF <- survfit(form, riskDF)
  new("plasma",
      meanPredictions = meanPreds,
      riskDF = riskDF,
      fullModel = model1,
      riskModel = rmodel,
      splitModel = smodel,
      SF = SF)
}

## summary method for plasma objects
setMethod("summary", "plasma", function(object, ...) {
})

## plot method for plasma objects
setMethod("plot", c("plasma", "missing"), function(x, y,  col = c("blue", "red"), lwd = 2, xlab = "", ylab = "Fraction Surviving", mark.time = TRUE, legloc = "topright", ...) {
  plot(x@SF, col = col, lwd = lwd, xlab = xlab, ylab = ylab, mark.time = mark.time, ...)
  legend(legloc, paste(c("low", "high"), "risk"), col = col, lwd = lwd)
})

## predict method for plasma objects
setMethod("predict", "plasma", function(object, newdata = NULL,
                                             type = c("components", "risk", "split"),
                                             ...) {
  tempor <- lapply(names(trainData), function(N) {
    cat(N, "\n", file = stderr())
    wb <- wickedBig[[N]]
    goForIt <- lapply(names(trainData), function(M) {
      cat("\t", M, "\n", file = stderr())
      localModel <- wickedBig[[N]]$mango[[M]]$learn
      localTest <- t(testData[[N]])[, dimnames(localModel$coefficients)[[1]]]
      if (!is.null(localModel)) {
        predictions <- predict(localModel, localTest)
      } else {
        predictions <- NULL
      }
      predictions
    })
    names(goForIt) <- names(trainData)
    ## reintegrate
    lap <- lapply(goForIt, function(G) {
      top <- dim(G)[3]
      G[,,top]
    })
    do.call(cbind, lap)
  })
  names(tempor) <- names(trainData)

  
  tstArray <- array(NA, dim = c(nrow(testOut), ncol(meanPreds), 9))
  dimnames(tstArray) <- list(rownames(testOut), colnames(meanPreds), names(trainData))
  for (I in 1:length(tempor)) {
    B <- tempor[[I]]
    w <- colnames(B) == "epitype"
    if (any(w)) colnames(B)[w] <- "epitype1"
    tstArray[rownames(B), colnames(B), I] <- B
  }
  meanTestPreds <- as.data.frame(apply(tstArray, 1:2, mean, na.rm = TRUE))
  dim(meanTestPreds)
  all(!is.na(meanTestPreds))
  
  
  testcox <- predict(mod1, meanTestPreds)
  testOut$testCox <- testcox
  testOut$Split <- factor(c("low", "high")[1+1*(testcox > median(barf))], levels = c("low", "high"))
  
  testMod <- coxph(Surv(survyr, sstat == "dead") ~ Split, data = testOut)
  testMod
  coxph(Surv(survyr, sstat == "dead") ~ testCox, data = testOut)
  S <- summary(testMod)
  PT <- S$sctest[3]
  text(15.5, 0.75, paste("p =", formatC(PT, format = "e", digits = 2)))
})
