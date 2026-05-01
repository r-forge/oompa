### (C) Copyright 2026 Polina Bombina and Kevin R. Coombes
###
### Functions to compute path metrics
### Simplified into its own "source-able" file by KRC

###################################
### Length Distortion
###
path_length <- function(coords) {
  if (nrow(coords) < 2) return(NA_real_)
  sum(sqrt(rowSums(diff(coords)^2)))
}
LengthDistortion <- function(M1, M2) {
  log(path_length(M1) / path_length(M2))
}

###################################
### Segment Variance
###
segment_lengths <- function(coords) {
  if (nrow(coords) < 2) return(numeric(0))
  var(sqrt(rowSums(diff(coords)^2)))
}
SegmentVariance <- function(M1, M2) {
  log(segment_lengths(M1) / segment_lengths(M2))
}

###################################
### Curvature
###
compute_curvature <- function(coords) {
  n <- nrow(coords)
  if (n < 3) return(rep(NA_real_, n))
  p1 <- coords[-c(n - 1, n), , drop = FALSE]
  p2 <- coords[-c(1, n),     , drop = FALSE]
  p3 <- coords[-c(1, 2),     , drop = FALSE]
  v1 <- p2 - p1
  v2 <- p3 - p2
  cosang <- rowSums(v1 * v2) / (sqrt(rowSums(v1^2)) * sqrt(rowSums(v2^2)))
  ang    <- acos(pmin(pmax(cosang, -1), 1))
### huh?  c(NA_real_, ang, NA_real_)
  return(mean(ang))
}
Curvature <- function(M1, M2) {
  log(compute_curvature(M1) / compute_curvature(M2))
}

###################################
### Spatial Similarity
###
compute_spatial_similarity <- function(coords) {
  if (nrow(coords) < 2) return(NA_real_)
  dist_mat <- as.matrix(dist(coords))
  sim_mat  <- 1 / (1 + dist_mat)
  mean(sim_mat[upper.tri(sim_mat)])
}
SpatialSimilarity <- function(M1, M2) {
  log(compute_spatial_similarity(M1) / compute_spatial_similarity(M2))
}

###################################
### KRC: do we need this?
###
compute_metrics_one_path <- function(path_obj) {
  Xp <- path_obj$hd_coords
  Yp <- path_obj$ld_coords
  
  L_hd           <- path_length(Xp)
  L_2d           <- path_length(Yp)
  distortion_log <- if (!is.na(L_hd) && L_hd > 0) log(L_2d / L_hd) else NA_real_
  
  curv_hd       <- compute_curvature(Xp)
  curv_2d       <- compute_curvature(Yp)
  mean_curv_hd  <- mean(curv_hd, na.rm = TRUE)
  mean_curv_2d  <- mean(curv_2d, na.rm = TRUE)
  curvature_log <- if (!is.na(mean_curv_hd) && mean_curv_hd > 0) log(mean_curv_2d / mean_curv_hd) else NA_real_
  
  seg_hd      <- segment_lengths(Xp)
  seg_2d      <- segment_lengths(Yp)
  var_hd      <- if (length(seg_hd) > 1) var(seg_hd) else NA_real_
  var_2d      <- if (length(seg_2d) > 1) var(seg_2d) else NA_real_
  seg_var_log <- if (!is.na(var_hd) && var_hd > 0) log(var_2d / var_hd) else NA_real_
  
  spatial_hd  <- compute_spatial_similarity(Xp)
  spatial_2d  <- compute_spatial_similarity(Yp)
  spatial_log <- if (!is.na(spatial_hd) && spatial_hd > 0) log(spatial_2d / spatial_hd) else NA_real_
  
  data.frame(
    n_clusters      = nrow(Yp),
    length_hd       = L_hd,
    length_2d       = L_2d,
    distortion_log  = distortion_log,
    mean_curv_hd    = mean_curv_hd,
    mean_curv_2d    = mean_curv_2d,
    curvature_log   = curvature_log,
    seg_var_hd      = var_hd,
    seg_var_2d      = var_2d,
    seg_var_log     = seg_var_log,
    spatial_sim_hd  = spatial_hd,
    spatial_sim_2d  = spatial_2d,
    spatial_sim_log = spatial_log,
    stringsAsFactors = FALSE
  )
}
