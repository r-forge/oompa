
setClass("CombinedWeights",
         slots = c(combined = "matrix",
                   featureSize = "numeric",
                   dataSource = "factor"))

combineAllWeights <- function(pl) {
  combined <- lapply(names(pl@compModels), function(N) getAllWeights(pl, N)@contrib)
  names(combined) <- names(pl@compModels)
  featsize <- sapply(combined, nrow)
  datasrc <- factor(rep(names(combined), times = featsize))
  combined <- do.call(rbind, combined)
  new("CombinedWeights",
      combined = combined,
      featureSize = featsize,
      dataSource = datasrc)
}

## summary method for CombinedWeights objects
setMethod("summary", "CombinedWeights", function(object, ...) {
  contra <- object@combined
  temp <- list(object@dataSource)
  aggie <- function(X, L, F) {
    mu <- aggregate(X, L, F)
    rownames(mu) <- mu[, 1]
    mu$Group.1 <- NULL
    mu <- as.matrix(mu)
  }

  mean <- aggie(contra, temp, mean)
  sd <- aggie(contra, temp, sd)
  med <- aggie(contra, temp, median)
  mad <- aggie(contra, temp, mad)

  list(Mean = mean, SD = sd, Median = med, MAD = mad)
})
