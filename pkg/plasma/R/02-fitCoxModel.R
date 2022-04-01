setOldClass("plsRcoxmodel")
setOldClass("coxph")
setOldClass("survfit")

setClass("SingleModel",
         slots=c(plsmod = "plsRcoxmodel",
                 Xout = "data.frame",
                 SF = "survfit",
                 riskModel = "coxph",
                 splitModel = "coxph"))

fitSingleModel <- function(object, N, timevar, eventvar, eventvalue) {
  ## get the training data
  X <- object@data[[N]]
  allNA <- apply(X, 2, function(xcol) all(is.na(xcol)))
  X <- X[, !allNA]
  ident <- apply(X, 1, function(x) length(unique(x)))
  X <- X[ident > 1, ]
  out <-  object@outcome
  Xout <-out[colnames(X),]
  mynt <- round(1 + log10(nrow(X)))
  ## Fit the PLS model
  plsmod <- plsRcoxmodel(t(X), time = Xout[, timevar],
                         event = Xout[eventvar] == eventvalue,
                         nt = mynt)
  Xout$Risk = predict(plsmod)
  Xout$Split <- 1*(Xout$Risk > median(Xout$Risk))
  riskModel <- coxph(formula(paste("Surv(", timevar, ",", eventvar, "== \"",
                                   eventvalue, "\") ~ Risk", sep = "")),
                     data = Xout)
  splitModel <- coxph(formula(paste("Surv(", timevar, ",", eventvar, "==\"",
                                    eventvalue, "\") ~ Split", sep = "")),
                      data = Xout)
  SF <- survfit(formula(paste("Surv(", timevar, ",", eventvar, "==\"",
                                    eventvalue, "\") ~ Split", sep = "")),
                      data = Xout)
  new("SingleModel",
      plsmod = plsmod,
      Xout = Xout,
      SF = SF,
      riskModel = riskModel,
      splitModel = splitModel)
}


setMethod("predict", "SingleModel", function(object, newdata = NULL,
                                             type = c("components", "risk", "split"),
                                             ...) {
  type <- match.arg(type)
  model <- switch(type,
                  components = object@plsmod,
                  risk = object@riskModel,
                  split = object@splitModel)
  arglist = list(model)
  if(!is.null(newdata)) arglist <- c(arglist, newdata)
  do.call(predict, arglist)
})

## summary method for SingleModel objects
setMethod("summary", "SingleModel", function(object, ...) {
  cat("Risk Model:\n", file = stdout())
  print(summary(object@riskModel))
  cat("SplitModel:\n", file = stdout())
  print(summary(object@splitModel))
  cat("PLS Model:\n", file = stdout())
  S <- summary(object@plsmod$FinalModel)
  S$call <- NULL
  S
})

setMethod("plot", c("SingleModel", "missing"), function(x, y,  col = c("blue", "red"), lwd = 2, xlab = "", ylab = "Fraction Surviving", mark.time = TRUE, legloc = "topright", ...) {
  plot(x@SF, col = col, lwd = lwd, xlab = xlab, ylab = ylab,  mark.time = mark.time, ...)
  legend(legloc, paste(c("Low", "High"), "Risk"), col = col, lwd = lwd)
})



fitCoxModels <- function(object, timevar, eventvar, eventvalue, verbose = TRUE) {
  firstPass <- lapply(names(object@data), function(N) {
    if(verbose) cat(N, "\n", file = stderr())
    fitSingleModel(object, N)
  })
  names(firstPass) <- names(object@data)
  firstPass
}

getSizes <- function(object) {
  NT <- sapply(firstPass, function(result) {
  S <- summary(result$Qmodel)
  PT <- S$sctest[3]
  c(NT = result$plsmod$nt, cNT = result$plsmod$computed_nt, p = PT)
  })
  colnames(NT) <- names(firstPass)
  t(NT)
 }

showme <- function(object) {
  S <- summary(object$Qmodel)
  PT <- S$sctest[3]
  ## Plot the results
  opar <- par(mfrow = c(1, 2))
  plot(object$Xsf, 
       col=c("red", "blue"), lwd = 2,
       main = paste("Cox model on training data (", object$N, ")", sep = ""))
  legend("topright", paste(c("Low", "High"), "Risk"), col = c("red", "blue"), lwd=2)
  plot(object$Qsf, 
       col=c("red", "blue"), lwd = 2,
       main = paste("Cox model on test data (", object$N, ")", sep = ""))
  legend("topright", paste(c("Low", "High"), "Risk"), col = c("red", "blue"), lwd=2)
  text(15.5, 0.75, paste("p =", formatC(PT, format = "e", digits = 2)))
  par(opar)
}

showAll <- function(object) {
  sapply(firstPass, showme)
}

