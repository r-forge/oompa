###################################
### Average Jaccard Distance
###

# Function to calculate the Jaccard distance between two sets
jaccard_distance <- function(setA, setB) {
  intersection_size <- length(intersect(setA, setB))
  union_size <- length(union(setA, setB))
  jaccard_dist <- 1 - (intersection_size / union_size)
  return(jaccard_dist)
}

# Function to compute Average Jaccard Distance
AverageJaccardDistance <- function(M1, M2, k = 50) {
  ## Get the number of points (cells)
  num_points <- ncol(M1)
  ## Ensure k is not greater than the number of points
  k <- min(k, num_points - 1)
  ## Compute k-nearest neighbors in the original space
  original_neighbors <- get.knn(M1,  k = k + 1)$nn.index
  ## Compute k-nearest neighbors in the embedding space
  embedding_neighbors <- get.knn(M2,  k = k + 1)$nn.index
  ## Initialize vector to store Jaccard distances
  jaccard_distances <- numeric(num_points)
  ## Calculate Jaccard distance for each point
  for (i in 1:num_points) {
    setA <- original_neighbors[i, -1]  ## Exclude the point itself
    setB <- embedding_neighbors[i, -1]  ## Exclude the point itself
    jaccard_distances[i] <- jaccard_distance(setA, setB)
  }
  ## Compute the average Jaccard distance
  average_jaccard_distance <- mean(jaccard_distances, na.rm = TRUE)
  return(average_jaccard_distance)
}

###################################
### Continuity and Trustworthiness
###
compute_ct_list <- function(M1, M2, lastNeighbor = 300) {
  ## Perform ContTrustMeasure
  ct_pca_result <- ContTrustMeasure(datamat = M1,
                                    projmat = M2,
                                    lastNeighbor = lastNeighbor)
    ## Return the result
    return(ct_pca_result)
}
AvgContinuity <- function(M1, M2, k = 300) {
  ct <- compute_ct_list(M1, M2, lastNeighbor = k)
  return(ct[,6])
}
AvgTrustworthiness <- function(M1, M2, k = 300) {
  ct <- compute_ct_list(M1, M2, lastNeighbor = k)
  return(ct[,3])
}


###################################
### Co-Ranking Matrix
###
QNX <- function(D1, D2) {
  ## Convert to matrices
  D1_matrix <- .matrify(D1)
  D2_matrix <- .matrify(D2)
  ## Compute the co-ranking matrix and Q_NX
  q_value <- coranking(D1_matrix, D2_matrix)
  qnx <- Q_NX(q_value)
  return(qnx)
}
