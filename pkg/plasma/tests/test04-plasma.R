library(plasma)
data("TCGA-ESCA")
## prepare MultiOmics
MO <- prepareMultiOmics(assemble, Outcome)
## test complete cox models
bigfit <- fitCoxModels(MO, "Days", "vital_status", "dead")
## extend across dataset pairs
extension <- extendCoxModels(MO, bigfit)
class(extension)
dim(extension)
heatmap(extension, scale = "none")


mfm <- plasma(MO, bigfit)
