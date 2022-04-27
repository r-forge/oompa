library(plasma)
data("TCGA-ESCA")
#assemble <- assemble[-1] # remove clinical data until we convert to numeric
MO <- prepareMultiOmics(assemble, Outcome)
train <- rep(c(TRUE, FALSE), times = c(112, 185-112))
MO2 <- MO[, train]
## Chweck each dataset
if (FALSE) {
  fitted <- fitSingleModel(MO, "miRSeq", "Days", "vital_status", "dead")
  summary(fitted) # OK
  fitted <- fitSingleModel(MO, "Meth450", "Days", "vital_status", "dead")
  summary(fitted) # OK
  fitted <- fitSingleModel(MO, "mRNASeq", "Days", "vital_status", "dead")
  summary(fitted) # OK
  fitted <- fitSingleModel(MO, "RPPA", "Days", "vital_status", "dead")
  summary(fitted) # OK
}
fitted <- fitSingleModel(MO2, "MAF", "Days", "vital_status", "dead")
summary(fitted)

plot(fitted, xlab = "Time (Days)", legloc = "topright", main = "Training Data")
p <- predict(fitted)
p <- try( predict(fitted, "risk") )
p <- try( predict(fitted, type = "riak") )
p <- predict(fitted, type = "risk")
q <- predict(fitted, type = "split")
plot(p, q)

testobj <- MO[, !train]
summary(testobj)
pre <- predict(fitted, newdata = testobj)
prer <- predict(fitted, newdata = testobj, type = "risk")
pres <- predict(fitted, newdata = testobj, type = "split")
pairs(cbind(pre, prer, pres))
