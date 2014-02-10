library(ClassComparison)

set.seed(782562)
betaPerc <- 0.10
nGenes <- 1000
alpha <- 0.2
nUnif <- trunc(nGenes*(1-betaPerc))
pvals <- c(runif(nUnif), rbeta(nGenes-nUnif, alpha, 1))
fit <- Bum(pvals)

fit@ahat
fit@lhat
fit@pihat

fdr <- 0.2
countSignificant(fit, alpha=fdr, by="FDR")
picked <- selectSignificant(fit, alpha=fdr, by="FDR")

truth <- rep(c("Null", "Real"), times=c(nUnif, nGenes-nUnif))

table(truth, picked)
