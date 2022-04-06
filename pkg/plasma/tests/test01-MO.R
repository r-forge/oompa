library(plasma)
data("TCGA-ESCA")
ls()
MO <- prepareMultiOmics(assemble, Outcome)
summary(MO)
opar <- par(mai = c(1.02, 1.82, 0.82, 0.42))
plot(MO)
par(opar)

MO2 <- MO[, 1:50]
summary(MO2)

## disorient one of the component datasets.
badinput <- assemble
badinput$mRNASeq <- t(badinput$mRNASeq)
pmo <- try( prepareMultiOmics(badinput, Outcome) )
