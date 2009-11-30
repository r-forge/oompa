## bimIndex.R

bimodalIndex <- function(dataset, verbose=TRUE) {
  bim <- matrix(NA, nrow=nrow(dataset), ncol=6)
  if (verbose) cat("1 ")
  for (i in 1:nrow(dataset)) {
    if (verbose && 0 == i%%100) cat(".")
    if (verbose && 0 == i%%1000) cat(paste("\n", 1 + i/1000, ' ', sep=''))
    x <- as.vector(as.matrix(dataset[i, ]))
    if (any(is.na(x))) next
    mc <- Mclust(x, G = 2, modelNames = "E")
    sigma <- sqrt(mc$parameters$variance$sigmasq)
    delta <- abs(diff(mc$parameters$mean))/sigma
#    pi <- max(mc$parameters$pro)
    pi <- mc$parameters$pro[1]
    bi <- delta * sqrt(pi*(1-pi))
    bim[i,] <-  c(mc$parameters$mean, sigma=sigma, delta=delta, pi=pi, bim=bi)
  }
  if(verbose) cat("\n")
  bim <- as.data.frame(bim)
  dimnames(bim) <- list(rownames(dataset),
                        c("mu1", "mu2", "sigma", "delta", "pi", "BI"))
  bim
}
