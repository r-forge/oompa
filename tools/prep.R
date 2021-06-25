### sometimes emacs is stupid
plot(2,3)
dev.off()

options(install.packages.compile.from.source = "never")

#source("https://bioconductor.org/biocLite.R")
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install(update = FALSE, ask = FALSE)
if (!require("Biobase")) {
#  biocLite(suppressUpdates = TRUE, suppressAutoUpdate = TRUE)
  BiocManager::install("Biobase", update = FALSE, ask = FALSE)
}
mybiocs <- c("edgeR",
             "Biostrings",
             "affy",
             "geneplotter",
             "DNAcopy",
             "flowCore",
             "AnnotationDbi",
             "RSQLite"
             )
for (p in mybiocs) {
  if (!require(p, character.only=TRUE)) {
    BiocManager::install(p, update = FALSE, ask = FALSE)
#    biocLite(p, suppressUpdates = TRUE, suppressAutoUpdate = TRUE)
  }
}

mypacks <- c("fortunes",
             "knitr",
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
             "XML", # not available for 3.6.3; copied from 3.6.0?
             "nFactors",
             "NbClust",
             "fgui",
             "epiR",
             "scatterplot3d",
             "quantmod",
             "gtools",
             "ggplot2",
             "scales",
             "plyr",
             "sirt",
             "Rtsne",
             "igraph",
             "dendextend",
             "DirichletReg",
             "rjson",
             "kohonen",
             "umap",
             "isotone",
             "fields",
             "flexmix" # why?
             )
for (p in mypacks) {
  if (!require(p, character.only = TRUE, quietly = TRUE)) {
    install.packages(p)
    library(p, character.only=TRUE)
  }
}
