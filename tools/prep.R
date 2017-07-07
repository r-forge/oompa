### sometimes emacs is stupid
plot(2,3)

source("https://bioconductor.org/biocLite.R")
if (!require("Biobase")) {
  biocLite()
}
mybiocs <- c("edgeR",
             "Biostrings",
             "affy")
for (p in mybiocs) {
  if (!require(p, character.only=TRUE)) {
    biocLite(p)
  }
}

mypacks <- c("knitr",
             "rmarkdown",
             "RColorBrewer",
             "colorspace",
             "rgl",
             "vioplot",
             "xtable",
             "mclust",
             "e1071",
             "randomForest",
             "neuralnet",
             "mc2d",
             "combinat",
             "abind",
             "ltm",
             "foreach",
             "doParallel",
             "KernSmooth",
             "matrixStats",
             "kernlab",
             "changepoint",
             "cpm",
             "ade4",
             "mgcv",
             "quantreg",
             "robustbase",
             "cobs",
             "timeDate",
             "drc",
             "DoseFinding",
             "alabama",
             "doSNOW",
             "movMF",
             "XML"
             )
for (p in mypacks) {
  if (!require(p, character.only=TRUE)) {
    install.packages(p)
    library(p, character.only=TRUE)
  }
}
