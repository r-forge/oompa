setClass("plastron",
         slots = c(train = "MultiOmics",
                   test = "MultiOmics",
                   timeVariable = "character",
                   eventVariable = "character"))

doIt <- function(N, object) {
  ## get the training data
  X <- trainData[[N]]
  allNA <- apply(X, 2, function(xcol) all(is.na(xcol)))
  X <- X[, !allNA]
  ident <- apply(X, 1, function(x) length(unique(x)))
  X <- X[ident > 1, ]
  Xout <- trainOut[colnames(X),]
  mynt <- round(1 + log10(nrow(X)))
  ## get the test data
  Q <- testData[[N]]
  allNA <- apply(Q, 2, function(xcol) all(is.na(xcol)))
  Q <- Q[, !allNA]
  Q <- Q[ident > 1,]
  Qout <- testOut[colnames(Q),]
  ## check trainng predictions
  opera <- predict(coxmod$FinalModel)
  Xout$Split <- factor(c("Low", "High")[1 + 1*(opera > median(opera))],
                       levels = c("Low", "High"))
  Xmodel <- coxph(Surv(survyr, sstat == "dead") ~ Split, data = Xout)
  Xsf <- survfit(Surv(survyr, sstat == "dead") ~ Split, data = Xout)
  ## Check test presdictions
  opera2 <- predict(coxmod, t(Q))
  Qout$Split <- opera2 > median(opera) # use training split point
  Qmodel <- coxph(Surv(survyr, sstat == "dead") ~ Split, data = Qout)
  Qsf <- survfit(Surv(survyr, sstat == "dead") ~ Split, data = Qout)

  ## save something
  invisible(list(N=N, plsmod = coxmod, Xmodel = Xmodel, Xsf = Xsf, Qmodel = Qmodel, Qsf = Qsf))
}

setOldClass("plsRcoxmodel")
setOldClass("coxph")
setOldClass("survfit")

setClass("SingleModel",
         slots=c(plsmod = "plsRcoxmodel",
                 Xout = "data.frame",
                 SF = "survfit",
                 riskModel = "coxph",
                 splitModel = "coxph"))

fitOneCoxModel <- function(object, N, timevar, eventvar, eventvalue) {
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

predictOneCoxModel <- function(model, newdata) {
}



fitCoxModels <- function(object, timevar, eventvar, eventvalue) {
  firstPass <- lapply(willuse, function(N) {
    cat(N, "\n", file = stderr())
    doIt(N)
  })
  names(firstPass) <- willuse
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

