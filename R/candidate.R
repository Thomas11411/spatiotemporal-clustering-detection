#' Build Candidate Spatial Clusters
#'
#' @description Generates all possible candidate spatial cluster combinations based on 
#' given geographic coordinates and a maximum distance threshold, and transforms them into 
#' an indicator matrix (Lambda) for subsequent spatiotemporal scan and cluster detection.
#'
#' @param location A matrix or data frame containing the geographic coordinates of each spatial point (e.g., latitude/longitude or 2D planar coordinates).
#' @param max.dist Numeric value representing the maximum allowable distance threshold for forming a cluster.
#'
#' @return A sparse or binary indicator matrix \code{Lambda} of dimension \code{n x num_unique}, 
#' where rows represent spatial points and columns represent unique candidate cluster combinations 
#' (1 indicates the point belongs to the cluster, 0 otherwise).
#' 
#' @export
candidate <- function(location, max.dist) {
  # Calculate the Euclidean distance matrix between all spatial points
  dist_mat <- as.matrix(dist(location))
  n <- nrow(location)
  
  cluster_list <- vector("list", n)
  
  # For each spatial point, find neighboring points within max.dist and order them by distance to build progressive clusters
  for (i in 1:n) {
    dists <- dist_mat[i, ]
    in_range <- which(dists < max.dist)
    ord <- in_range[order(dists[in_range])] # Sort only points within max.dist
    
    # Build progressive cluster index lists for center point i (e.g., expanding from 1 point to multiple)
    cluster_list[[i]] <- lapply(seq_along(ord), function(j) ord[1:j])
  }
  
  # Unlist and remove duplicate cluster combinations (sort node sets before global deduplication)
  all_clusters <- unlist(cluster_list, recursive = FALSE)
  sorted_clusters <- lapply(all_clusters, sort)
  unique_clusters <- unique(sorted_clusters)
  
  num_unique <- length(unique_clusters)
  Lambda <- matrix(0, nrow = n, ncol = num_unique)
  
  # Populate the indicator matrix Lambda with unique cluster combinations
  for (k in seq_len(num_unique)) {
    Lambda[unique_clusters[[k]], k] <- 1
  }
  
  return(Lambda)
}
