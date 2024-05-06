### sometimes emacs is stupid
plot(2,3)
dev.off()

options(install.packages.compile.from.source = "never")

if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install(update = FALSE, ask = FALSE)
if (!require("Biobase")) {
  BiocManager::install("Biobase", update = FALSE, ask = FALSE)
}
mybiocs <- c("edgeR",
             "Biostrings",
             "affy",
             "geneplotter",
             "DNAcopy",
             "flowCore",
             "AnnotationDbi",
             "RSQLite",
             "mixOmics",
             "msa",
             "survcomp",
             "DescTools",
             "org.Hs.eg.db",
             "graph",
             "Rgraphviz" 
             )
for (p in mybiocs) {
  if (!require(p, character.only=TRUE)) {
    BiocManager::install(p, update = FALSE, ask = FALSE)
  }
}

mypacks <- c("fortunes",
             "V8",
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
             "rjson",
             "kohonen",
             "umap",
             "isotone",
             "fields",
             "flexmix", # why?
             "pls",
             "plsRcox",
             "ape",
             "corrplot",
             "DirichletReg",
             "R.rsp",
             "beanplot",
             "circlize",
             "TDA",
             "PubChemR"
             )
for (p in mypacks) {
  if (!require(p, character.only = TRUE, quietly = TRUE)) {
    install.packages(p)
    library(p, character.only=TRUE)
  }
}
