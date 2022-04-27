library(plasma)
data("TCGA-ESCA")
ls()
MO <- prepareMultiOmics(assemble, Outcome)
summary(MO)
opar <- par(mai = c(1.02, 1.82, 0.82, 0.42))
plot(MO)
par(opar)

train <- rep(c(TRUE, FALSE), times = c(112, 185-112))
MO2 <- MO[, train]
summary(MO2)

if (require("survival")) {
  plot(survfit(Surv(Days, vital_status == "dead") ~ train, data = MO@outcome))
}


## disorient one of the component datasets.
badinput <- assemble
badinput$mRNASeq <- t(badinput$mRNASeq)
pmo <- try( prepareMultiOmics(badinput, Outcome) )
