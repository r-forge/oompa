library(plasma)
data("TCGA-ESCA")
## prepare MultiOmics
MO <- prepareMultiOmics(assemble, Outcome)
## test complete cox models
set.seed(12345)
train <- rep(FALSE, 185)
train[sample(185, 113)] <- TRUE
MO2 <- MO[, train]
bigfit <- fitCoxModels(MO2, timevar = "Days",
                       eventvar = "vital_status",
                       eventvalue = "dead")
class(bigfit)
summary(bigfit)
getSizes(bigfit)
plot(bigfit)
preds <- predict(bigfit)
