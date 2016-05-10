source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("Biostrings")
biocLite("edgeR")

deps <- c("kernlab", "mclust", "cobs", "rpart", "ltm", "doParallel",
          "colorspace", "timeDate", "quantreg", "nnet", "e1071",
          "cluster", "abind", "foreach", "rgl", "RColorBrewer",
          "xtable", "mc2d", "mgcv", "class", "randomForest",
          "combinat", "KernSmooth", "ade4", "changepoint", "movMF",
          "robustbase", "cpm")

install.packages(deps)

if (FALSE) {
  for (x in c("Biobase", "Biostrings", deps)) {
    library(x, character.only=TRUE)
  }
}
