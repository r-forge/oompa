#library(FNN)  # For k-nearest neighbors

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

# Function to compute Average Jaccard Distance for a specific repetition and method
compute_average_jaccard_distance <- function(original_data, embedding_data, k) {
  ## Get the number of points (cells)
  num_points <- ncol(original_data)
  ## Ensure k is not greater than the number of points
  k <- min(k, num_points - 1)
  ## Compute k-nearest neighbors in the original space
  original_neighbors <- get.knn(original_data,  k = k + 1)$nn.index
  ## Compute k-nearest neighbors in the embedding space
  embedding_neighbors <- get.knn(embedding_data,  k = k + 1)$nn.index
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
