##

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

setMethod("image", "CombinedWeights", function(x, ...) {
  summer <- summary(x)
  mx <- max(abs(summer$Mean))
  sx <- max(summer$SD)
  dx <- max(abs(summer$Median))
  ax <- max(summer$MAD)
  N <- ncol(summer$Mean)
  M <- nrow(summer$Mean)
  opar <- par(mfrow = c(2,2))
  on.exit(par(opar))
  image(1:N, 1:M, t(summer$Mean), xlab = "Component", ylab = "Dataset",
        main = "Mean", zlim = c(-mx, mx), col = redgreen(64))
  image(1:N, 1:M, t(summer$SD), xlab = "Component", ylab = "Dataset",
        main = "SD", zlim = c(0, sx), col = bluescale(64))
  image(1:N, 1:M, t(summer$Median), xlab = "Component", ylab = "Dataset",
        main = "Median", zlim =c(-dx, dx), col = redgreen(64))
  image(1:N, 1:M, t(summer$MAD), xlab = "Component", ylab = "Dataset",
        main = "Mad", zlim = c(0, ax),col = bluescale(64))
  invisible(x)
})

stdize <- function(object, type = c("standard", "robust")) {
  type <- match.arg(type)
  summer <- summary(object)
  mu <- switch(type,
               standard = summer$Mean,
               robust = summer$Median)
  sigma <- switch(type,
                  standard = summer$SD,
                  robust = summer$MAD)
  M <- apply(mu, 2, function(X) rep(X, times = object@featureSize))
  S <- apply(sigma, 2, function(X) rep(X, times = object@featureSize))
  (object@combined - M)/S
}

interpret <- function(object, component, alpha = 0.05) {
  brute <- stdize(object)
  Q <- qnorm(1 - alpha/2) # two-sided 5% cutoff
  topA <- which(abs(brute[, component]) > Q)
  data.frame(Feature = rownames(brute)[topA],
             Source = object@dataSource[topA],
             Weight = brute[topA, component])
}

# merge feature contributiosn across all components
getFinalWeights <- function(object) {
  fm <- object@fullModel
  coy <- matrix(fm$coefficients, ncol = 1)
  mainterms <- attr(terms(fm), "term.labels")
  combined <- lapply(names(object@compModels), function(N) {
    getAllWeights(object, N)@contrib[, mainterms]
  })
  names(combined) <- names(object@compModels)
  runthrough <- lapply(names(combined), function(N) {
    X <- combined[[N]] %*% coy
    data.frame(Weight = X, Source = N, Feature = rownames(X))
  })
  FW <- do.call(rbind, runthrough)
  mu <- aggregate(FW$Weight, list(FW$Source), mean)
  mu2 <- rep(mu$x, times = as.vector(table(FW$Source)))
  sigma <- aggregate(FW$Weight, list(FW$Source), sd)
  sigma2 <- rep(sigma$x, times = as.vector(table(FW$Source)))
  FW$Standard <- (FW$Weight - mu2)/sigma2
  FW
}

setMethod("barplot", c("plasma"), function(height, source, n, ...) {
  wws <- getFinalWeights(height) # default name from function, but a plasma object
  wmut <- wws[wws$Source  == source,]
  ## Get positive and negative weights 
  wmut <- wws[wws$Source  == source, ]
  wmut <- wmut[order(wmut$Weight),]
  wmutt <- wmut[1:n,]
  wmutt <- wmutt[order(wmutt$Weight, decreasing = TRUE),]
  wmutt$Feature <- as.character(wmutt$Feature)
  wmutt$Feature <- factor(wmutt$Feature, levels = unique(wmutt$Feature))

  wmut2 <- wws[wws$Source  == source,]

  ## positive weights 
  wmut2 <- wmut2[order(wmut2$Weight, decreasing = TRUE),]
  wmutt2 <- wmut2[1:20,]
  wmutt2 <- wmutt2[order(wmutt2$Weight),]
  wmutt2$Feature <- as.character(wmutt2$Feature)
  wmutt2$Feature <- factor(wmutt2$Feature, levels = unique(wmutt2$Feature))
  wmutt <- rbind(wmutt, wmutt2)

  ## Define variable mapping
  map <- aes_string(x = "Weight", y = "Feature", fill = "Weight")

  ## Make the ggplot
  p <- ggplot(wmutt, map) + geom_bar(stat = "identity") +  theme_bw() + 
    scale_fill_continuous(low = "red", high = "blue") + 
    guides(fill = guide_colorbar(title = source , reverse = TRUE))
  ## return it
  p
})
