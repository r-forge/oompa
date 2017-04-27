source("https://bioconductor.org/biocLite.R")
if (!require("Biobase")) {
  biocLite()
}
mybiocs <- c("edgeR",
             "Biostrings")
for (p in mybiocs) {
  if (!require(p, character.only=TRUE)) {
    biocLite(p)
  }
}

### sometimrs emacs is stupid
plot(2,3)

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
             "movMF",
             "ade4",
             "mgcv",
             "quantreg",
             "robustbase",
             "cobs"
             )
for (p in mypacks) {
  if (!require(p, character.only=TRUE)) {
    install.packages(p)
  }
}
