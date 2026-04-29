###################################
### Milnor Distortion
###
setClass("Milnor",
         slot = c(scale = "numeric",
                  distortion = "numeric"))
Milnor <- function(D1, D2) {
  scale <- log(D2/D1)
  scale <- scale[!is.infinite(scale)]
  distortion <- diff(range(scale))
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

###################################
### Sigma Distortion
###
calculate_normalized_ratios <- function(original_dist_matrix, embedding_dist_matrix) {
  ## Avoid division by zero
  original_dist_matrix[original_dist_matrix == 0] <- 1e-10
  ## Compute ratio of distances
  rho_f <- embedding_dist_matrix / original_dist_matrix
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
compute_sigma_distortion <- function(original_dist_matrix, embedding_dist_matrix) {
  ## Compute normalized ratios
  rho_tilde_f <- calculate_normalized_ratios(original_dist_matrix, embedding_dist_matrix)
  ## Compute sigma-distortion
  sigma_distortion <- mean((rho_tilde_f - 1)^2)
  return(sigma_distortion)
}

###################################
### Stress
###
compute_stress <- function(original_dist_matrix, embedding_dist_matrix) {
  ## Calculate the squared differences between distances
  dist_diff_squared <- (original_dist_matrix - embedding_dist_matrix)^2
  ## Compute the sum of squared differences
  sum_dist_diff_squared <- sum(dist_diff_squared, na.rm = TRUE)
  ## Compute the sum of squared original distances
  sum_original_dist_squared <- sum(original_dist_matrix^2, na.rm = TRUE)
  ## Calculate stress
  stress <- sqrt(sum_dist_diff_squared / sum_original_dist_squared)
  return(stress)
}

###################################
### M1 Distortion (mean-squared)
###
frobenius_norm_squared <- function(matrix) {
  norm_squared <- norm(matrix, type = "F")^2
  return(norm_squared)
}

# Function to compute M1 distortion
compute_m1_distortion <- function(original_data, embedding_data) {
  ## Calculate Frobenius norms squared
  original_norm_squared <- frobenius_norm_squared(original_data)
  embedding_norm_squared <- frobenius_norm_squared(embedding_data)
  ## Calculate M1 distortion
  m1_distortion <- abs((embedding_norm_squared / original_norm_squared) - 1)
  return(m1_distortion)
}

###################################
### Spearman's Rho
###
### Uses the DRquality package
###
calc_spearman_rho <- function(original_dist, embedding_dist) {
  original_dist_matrix <- as.matrix(original_dist)
  method_dist_matrix <- as.matrix(method_dist)
  ## Calculate Spearman's Rho
  spearman_rho <- SpearmansRho(original_dist_matrix, method_dist_matrix)
  return(spearman_rho)
}


###################################
### Earth Mover's Distance
###
### Uses the emdist package
###
calc_emd <- function(original_dist, embedding_dist) {
  ## Convert to matrices
  original_dist_matrix <- as.matrix(original_dist)
  method_dist_matrix <- as.matrix(method_dist)
  ## Compute Earth Mover's Distance
  emdist_value <- emd(original_dist_matrix, method_dist_matrix, dist = "euclidean")
  return(emdist_value)
}
