library(plasma)
data("TCGA-ESCA")
MO <- prepareMultiOmics(assemble, Outcome)
fitted <- fitOneCoxModel(MO, "miRSeq", "Days", "vital_status", "dead")
## The next part should probably be made into 'summary' and 'plot' methods
## for the 'SingleModel' class

class(fitted)
slotNames(fitted)
summary(fitted@riskModel)
summary(fitted@splitModel)
summary(fitted@Xout)
plot(fitted@SF, col = c("blue", "red"), lwd = 2, xlab = "Days", mark.time = TRUE)
legend("bottomleft", paste(c("Low", "High"), "Risk"), col = c("blue", "red"),
       lwd=2)

