### (C) Copyright 2026 Polina Bombina and Kevin R. Coombes
###
### Functions to compute path metrics
### Simplified into its own "source-able" file by KRC

###################################
### Smoothness
###
comp_smooth <- function( coords) {
  ticks <- 1:nrow(coords) # arbitrary parametrization
  ## fit smooth curve to each coordinate as a function of the parameter
  fitted <- apply(coords, 2, function(X) predict(loess(X ~ ticks, span = 0.1)))
  delta <- fitted - coords # difference between actual path and smooth fitted curve
  SE <- apply(delta^2, 1, sum) # Euclidean distance from path to fitted curve
  sqrt(SE)
}
Smoothness <- function(M1, M2) {
  log(mean(comp_smooth(M1)) / mean(comp_smooth(M2)))
}
RankSmoothness <- function(M1, M2) {
  cor(comp_smooth(M1), comp_smooth(M2), method = "spearman")
}

###################################
### Coil
###
comp_coil <- function(coords) {
  troid <- apply(coords, 2, mean)
  temp <- rbind(coords, troid)
  ed <- dist(temp)
  L <- nrow(temp)
  dd <- as.matrix(ed)[L,][-L]
  diff(range(dd))
}
Coil <- function(M1, M2) {
  log(comp_coil(M1) / comp_coil(M2))
}

###################################
### Contact
###
comp_contact <- function(coords, q, type = c("mean", "max", "count")) {
  type <- match.arg(type)
  ed <- as.matrix(dist(coords)) # small values here mean close contact
  q10 <- quantile(ed, q)        # use a quantile to decide what small means
  ai <- which(ed < q10, arr.ind = TRUE) # get matrix coordinates of close cells
  spread <- abs(apply(ai, 1, diff)) # number of steps apart in the sequence
  val <- switch(type,
                mean = mean(spread),
                max = max(spread),
                count = sum(spread > 0.12*nrow(coords))) # MAGIC
  val
}
Contact <- function(M1, M2, q = 0.1, type = c("mean", "max", "count")) {
  type <- match.arg(type)
  log(comp_contact(M1, q, type) / comp_contact(M2, q, type))
}


###################################
### EndpointDisplacement
###
endpt_disp <- function(coords, k) {
  start <- apply(coords[1:k,], 2, mean)
  L <- nrow(coords)
  end <- apply(coords[(L-k+1):L, ], 2, mean)
  sqrt(sum((start-end)^2))
}
EndpointDisplacement <- function(M1, M2, k=5) {
  log(endpt_disp(M1, k) / endpt_disp(M2, k))
}

###################################
### CrossingNumber
###
segments_intersect <- function(p1, p2, p3, p4, eps = 1e-12) {
  orient <- function(a, b, c) {
    (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1])
  }
  o1 <- orient(p1, p2, p3)
  o2 <- orient(p1, p2, p4)
  o3 <- orient(p3, p4, p1)
  o4 <- orient(p3, p4, p2)
  (o1 * o2 < -eps) && (o3 * o4 < -eps)
}

crossing_number_raw <- function(coords) {
  n <- nrow(coords)
  if (n < 4) return(0)
  count <- 0
  for (i in 1:(n - 3)) {
    for (j in (i + 2):(n - 1)) {
      if (j == i + 1) next
      if (segments_intersect(coords[i,], coords[i+1,],
                             coords[j,], coords[j+1,])) {
        count <- count + 1
      }
    }
  }
  count
}

CrossingNumber <- function(M2, eps = 1e-8) {
  cn <- crossing_number_raw(M2)
  n <- nrow(M2)
  cn / ((n-2)*(n-1)/2 + eps)
}
