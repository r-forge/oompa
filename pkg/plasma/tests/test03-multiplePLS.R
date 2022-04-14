library(plasma)
data("TCGA-ESCA")
## prepare MultiOmics
MO <- prepareMultiOmics(assemble, Outcome)
## test complete cox models
bigfit <- fitCoxModels(MO, "Days", "vital_status", "dead")
class(bigfit)
summary(bigfit)
getSizes(bigfit)
plot(bigfit)
preds <- predict(bigfit)
