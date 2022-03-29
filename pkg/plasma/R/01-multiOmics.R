## Need a class to hold the training and test data

## A multiomics data set. The data is a list, with each entry a data set
## from a different assay. Like MOFA, should have the same samples
setClass("MultiOmics",
         slots = c(data = "list",
                   outcome = "data.frame"))
## As in MOFA, each data set must contain the same set of samples
validMultiOmics <- function(object) {
  sampleSizes <- sapply(object@data, ncol)
  okSS <- all(sampleSizes == nrow(object@outcome))
  if (okSS) {
    namesOK <- sapply(object@data, function(DS) {
      all(colnames(DS) == rownames(object@data))
    })
    valid <- ifelse(namesOK, TRUE, "Sample names in all datasets must agree.")
  } else {
    valid <- "All datasets must have the same number of samples."
  }
  valid
}
setValidity("MultiOmics", validMultiOmics)

## Here we assume that we have a list of datasets (on patient subsets)
## and a data.frame of outcomes thatn includes all patients.
## We assume that outcome is organized as patients x variables, but
## each dataset is organized as features x patients.
prepareMultiOmics <- function(datalist, outcome) {
  ready <- lapply(datalist, function(DS) {
    ## merge each dataset with the first two columns 
    filled <- merge(outcome[, 1:2], t(DS), by = "row.names", all.x = TRUE)
    ## first column should be the row names (patients)
    rownames(filled) <- filled[, 1]
    ## columns two and three are copies of the outcome data
    filled <- filled[, -(1:3)]
    t(filled)
  })
  names(ready) <- names(datalist)
  new("MultiOmics", data = ready, outcome = outcome)
}


