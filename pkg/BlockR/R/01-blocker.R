blockerLM <- function(dataset, sampleClades, proteinClades, NP, NS) {
  storage <- matrix(NA, nrow=NP*NS, ncol=3)
  idx <- 0
  for (np in 1:NP) {
    pgroup <- paste("P", cutree(proteinClades, k=np), sep='')
    for (ns in 1:NS) {
      sgroup <- paste("S", cutree(sampleClades, k=ns), sep='')
      pivot <- data.frame(Y=as.vector(as.matrix(dataset)),
                          P=rep(pgroup, each=nrow(dataset)),
                          S=rep(sgroup, times=ncol(dataset)))
      if (np == 1) {
        if (ns == 1) {
          model <- lm(Y ~ 1, data=pivot)
        } else {
          model <- lm(Y ~ S, data=pivot)
        }
      } else if (ns == 1) {
        model <- lm(Y ~ P, data=pivot)
      } else {
        model <- lm(Y ~ P*S, data=pivot)
      }
      idx <- idx+1
      storage[idx,] <- c(np, ns, mean(model$res^2)) # mean square residual
    }
  }
  colnames(storage) <- c("ProteinGroups", "SampleGroups", "MeanMSE")
  storage <- as.data.frame(storage)
  storage
}

## protein clades are columns
## sample clades are rows
blockerLoop <- function(dataset, sampleClades, proteinClades, NP, NS) {
  storage <- matrix(NA, nrow=NP*NS, ncol=3)
  idx <- 0
  for (np in 1:NP) {
    cat(np, "\n", file = stderr())
    for (ns in 1:NS) {
      RR <- factor(paste("R", cutree(sampleClades, k = ns), sep = ""))
      CC <- factor(paste("C", cutree(proteinClades, k = np), sep = ""))
      P <- 0*dataset
      for (R in levels(RR)) {
        for (C in levels(CC)) {
          P[RR == R, CC == C] <- mean(dataset[RR==R, CC==C])
        }
      }
      Res <- matrix(dataset - P, nrow = nrow(dataset))
      idx <- idx+1
      storage[idx,] <- c(np, ns, mean(Res^2)) # mean square residual
    } # end cn
  } # end rn
  colnames(storage) <- c("ProteinGroups", "SampleGroups", "MeanMSE")
  storage <- as.data.frame(storage)
  storage
}


blocker <- function(dataset) {
  sampleClades <- hclust(distanceMatrix(t(dataset), metric="pearson"), 
                         method="ward.D2")
  proteinClades <- hclust(distanceMatrix(dataset, metric="pear"), 
                          method="ward.D2")
  NP <- trunc(3 + sqrt(ncol(dataset)))
  NS <- trunc(1 + sqrt(nrow(dataset)))
  blockerLoop(dataset, sampleClades, proteinClades, NP, NS)
}

pivot <- function(results) {
  mat <- matrix(NA, nrow=max(results[,1]), ncol=max(results[,2]))
  for (i in 1:nrow(results)) {
    mat[results[i,1], results[i,2]] <- results[i,3]
  }
  mat
}

n2loglik <- function(x, nn) {
  xp <- x$ProteinGroups
  xs <- x$SampleGroups
  xe <- x$MeanMSE        # mean square residual
  kk <- xp*xs            # number of parameters
#  nn <- prod(dim(dat))   # number of observations
  n2loglik <- nn*log(xe) # -2 log(likelihood)
  aic <- n2loglik + 2*kk
  bic <- n2loglik + kk*log(nn)
  list(K=kk, N=nn, AIC=aic, BIC=bic, neg2ll=n2loglik)
}
