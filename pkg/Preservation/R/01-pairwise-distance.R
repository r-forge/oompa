##################################
### Utilities
###
.matrify <- function(D0) {
  if (inherits(D0, "dist")) {
    D0 <- as.matrix(D0)
  }
  D0
}


##################################
### Milnor Distortion
###
setClass("Milnor",
         slot = c(scale = "numeric",
                  distortion = "numeric"))
Milnor <- function(D1, D2) {
  D1 <- .matrify(D1)
  D2 <- .matrify(D2)
  scale <- log(D2/D1)
  scale <- scale[!is.infinite(scale)]
  distortion <- diff(range(scale, na.rm = TRUE))
  new("Milnor", scale = scale, distortion = distortion)
}
setMethod("summary", "Milnor", function(object, ...) {
  c(distortion = object@distortion,
    mean = mean(object@scale),
    median = median(object@scale),
    var = var(object@scale),
    iqr = IQR(object@scale))
})

if (!isGeneric("hist"))
  setGeneric("hist",
             function(x, ...) { standardGeneric("hist") }
             )

setMethod(hist, "Milnor",  function(x, main = "", xlab = "log(scale)", ...) {
  hist(x@scale, main = main, xlab = xlab, ...)
})

MilnorDistortion <- function(D1, D2) {
  M <- Milnor(D1, D2)
  M@distortion
}

###################################
### Sigma Distortion
###
calculate_normalized_ratios <- function(D1, D2) {
  ## Avoid division by zero
  D1[D1 == 0] <- 1e-10
  ## Compute ratio of distances
  rho_f <- D2 / D1
  ## Extract upper triangular part of the matrix
  rho_f <- rho_f[upper.tri(rho_f)]
  ## Compute normalization factor
  num_pairs <- length(rho_f)
  sum_rho_f <- sum(rho_f)
  ## Normalize distances
  rho_tilde_f <- (num_pairs * rho_f^2) / sum_rho_f
  return(rho_tilde_f)
}

# Function to compute sigma-distortion
SigmaDistortion <- function(D1, D2) {
  D1 <- .matrify(D1)
  D2 <- .matrify(D2)
  ## Compute normalized ratios
  rho_tilde_f <- calculate_normalized_ratios(D1, D2)
  ## Compute sigma-distortion
  sigma_distortion <- mean((rho_tilde_f - 1)^2)
  return(sigma_distortion)
}

###################################
### Stress
###
Stress <- function(D1, D2) {
  D1 <- .matrify(D1)
  D2 <- .matrify(D2)
  ## Calculate the squared differences between distances
  dist_diff_squared <- (D1 - D2)^2
  ## Compute the sum of squared differences
  sum_dist_diff_squared <- sum(dist_diff_squared, na.rm = TRUE)
  ## Compute the sum of squared original distances
  sum_D1_squared <- sum(D1^2, na.rm = TRUE)
  ## Calculate stress
  stress <- sqrt(sum_dist_diff_squared / sum_D1_squared)
  return(stress)
}

###################################
### M1 Distortion (mean-squared)
###
frobenius_norm_squared <- function(matrix) {
  norm_squared <- norm(matrix, type = "F")^2
  return(norm_squared)
}

## Function to compute M1 distortion
## Uses original data matrix, not dist
M1Distortion <- function(M1, M2) {
  ## Calculate Frobenius norms squared
  original_norm_squared <- frobenius_norm_squared(M1)
  embedding_norm_squared <- frobenius_norm_squared(M2)
  ## Calculate M1 distortion
  m1_distortion <- abs((embedding_norm_squared / original_norm_squared) - 1)
  return(m1_distortion)
}

###################################
### Spearman's Rho
###
### Uses the DRquality package
###
SpearmanRho <- function(D1, D2) {
  D1 <- .matrify(D1)
  D2 <- .matrify(D2)
  ## Calculate Spearman's Rho
  spearman_rho <- SpearmansRho(D1, D2)
  return(spearman_rho)
}

###################################
### Earth Mover's Distance
###
### Uses the emdist package
###
.process <- function(D1, res = 500) {
  R1 <- range(D1)
  delta1 <- diff(R1)/res/2
  breaks <- seq(min(D1) - delta1, max(D1) + delta1,length = res)
  H1 <- hist(D1, breaks = breaks, plot = FALSE)
  as.matrix(data.frame(H1$mids, H1$density))
}

EarthMover <- function(D1, D2) {
  ## Convert to distributions of distances
  H1 <- .process(D1)
  H2 <- .process(D2)
  ## Compute Earth Mover's Distance
  emdist_value <- emd(H1, H2, dist = "euclidean")
  return(emdist_value)
}

.target <- list(Stress = 0, MilnorDistortion = 0, M1Distortion = 0,
               SigmaDistortion = 0, SpearmanRho = 1, EarthMover = 0)
PWBest <- function(tag) {
  if (length(tag) > 1 ) {
    val <- unlist(.target[tag])
  } else {
    val <- .target[[tag]]
  }
  val
}
