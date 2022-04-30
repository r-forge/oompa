library(plasma)
data("TCGA-ESCA")
## prepare MultiOmics
MO <- prepareMultiOmics(assemble, Outcome)
MO <- MO[c("ClinicalBin", "ClinicalCont", "RPPA"),]
summary(MO)
## test complete cox models
bigfit <- fitCoxModels(MO, "Days", "vital_status", "dead")
## extend across dataset pairs
extension <- extendCoxModels(MO, bigfit)
class(extension)
dim(extension)
ext <- extension[!is.na(extension[,1]),]
heatmap(ext, scale = "none")


#mfm <- plasma(MO, bigfit)
