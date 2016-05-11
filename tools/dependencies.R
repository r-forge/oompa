source("https://bioconductor.org/biocLite.R")
biocLite(suppressUpdates=TRUE)
biocLite("Biostrings", suppressUpdates=TRUE)
biocLite("edgeR", suppressUpdates=TRUE)

deps <- c("kernlab", "mclust", "cobs", "rpart", "ltm", "doParallel",
          "colorspace", "timeDate", "quantreg", "nnet", "e1071",
          "cluster", "abind", "foreach", "rgl", "RColorBrewer",
          "xtable", "mc2d", "mgcv", "class", "randomForest",
          "combinat", "KernSmooth", "ade4", "changepoint", "movMF",
          "robustbase", "cpm")

install.packages(deps)

if (FALSE) {
  for (x in c("Biobase", "Biostrings", "edgeR", deps)) {
    library(x, character.only=TRUE)
    cat(x, "is available.\n", file=stderr())
  }
}
