#' Run Complete Spatiotemporal Cluster Detection and Analysis Pipeline
#'
#' @description Executes the full spatiotemporal analysis pipeline directly from an integrated 
#' shapefile containing time-series response variables, covariates, and coordinates. Automatically 
#' detects time steps from the shapefile attributes, and performs null model estimation, 
#' spatiotemporal correlation modeling, parallelized candidate cluster scanning, cluster 
#' significance testing, and iterative cluster regrouping.
#'
#' @param shp.file An \code{sf} spatial object (the integrated shapefile containing \code{Y_w*}, \code{X_w*}, \code{sp_X}, and \code{sp_Y} columns).
#' @param max.dist Numeric value, maximum distance threshold for building candidate clusters.
#' @param cutoff Numeric value, spatial distance cutoff for semivariogram estimation (defaults to 1.5).
#' @param T.cnt Numeric value, temporal trend degree (1 for intercept, 2 for linear, 3 for quadratic; defaults to 3).
#' @param nc Integer, number of CPU cores to use for parallel computing cluster (defaults to 11).
#'
#' @return A list containing the final fitted model object (\code{final.model}), 
#' the final regrouped cluster indicator matrix (\code{final.clusters}), 
#' the updated shapefile with cluster attributes (\code{updated.shp}), 
#' and the visualization plot object (\code{plot}).
#' 
#' @export
spatiotemporal_cluster_analysis <- function(shp.file, max.dist = 0.5, cutoff = 1.5, T.cnt = 3, nc = 11) {
  
  t0 <- Sys.time()
  
  # ------------------------------------------------------------------
  # Step 0: Automatically Extract Matrices and Coordinates from Shapefile
  # ------------------------------------------------------------------
  message("Step 0: Extracting response, covariate, and coordinates from integrated shapefile...")
  
  n_loc <- nrow(shp.file)
  
  # Automatically retrieve all column names starting with "Y_w" from the shapefile and sort them by time points
  y_cols <- sort(grep("^Y_w", names(shp.file), value = TRUE))
  Time <- length(y_cols)
  
  if (Time == 0) {
    stop("Error: No response variable columns starting with 'Y_w' were found in the shapefile attributes!")
  }
  
  # 1. Reconstruct response variable matrix
  response_mat <- matrix(NA, nrow = n_loc, ncol = Time)
  for (i in 1:Time) {
    response_mat[, i] <- shp.file[[y_cols[i]]]
  }
  
  # 2. Reconstruct covariate matrix (corresponding to X_w columns)
  x_cols <- sort(grep("^X_w", names(shp.file), value = TRUE))
  if (length(x_cols) != Time) {
    stop("Error: The number of covariate (X_w) columns does not match the number of response variable columns!")
  }
  
  covariate_mat <- matrix(NA, nrow = n_loc, ncol = Time)
  for (i in 1:Time) {
    covariate_mat[, i] <- shp.file[[x_cols[i]]]
  }
  
  # 3. Extract spatial coordinates
  coordinate <- data.frame(
    x = shp.file$sp_X / 1000,
    y = shp.file$sp_Y / 1000
  )
  
  # ------------------------------------------------------------------
  # Step 1: Generate Candidate Spatial Clusters
  # ------------------------------------------------------------------
  message("Step 1: Generating candidate spatial clusters...")
  Lambda.tmp <- candidate(coordinate, max.dist)
  
  response_critical <- colMeans(response_mat)
  group_counts <- t(Lambda.tmp == 1) %*% matrix(1, nrow = nrow(response_mat), ncol = 1)
  group_sums   <- t(Lambda.tmp == 1) %*% response_mat
  group_means  <- group_sums / as.vector(group_counts)
  
  threshold    <- round(Time / 2)
  condition    <- rowSums(sweep(group_means, 2, response_critical, ">")) > threshold
  i.set        <- which(condition)
  Lambda       <- Lambda.tmp[, i.set, drop = FALSE]
  nclust       <- ncol(Lambda)
  
  message(paste("-> Filtered candidate clusters count:", nclust))
  
  # ------------------------------------------------------------------
  # Step 2: Initialize Parallel Computing Cluster
  # ------------------------------------------------------------------
  message(paste("Step 2: Initializing parallel cluster with", nc, "cores..."))
  cl <- makeCluster(nc)
  
  clusterEvalQ(cl, {
    library("MASS")
    library("ggplot2")
    library("gstat")
    library("Matrix")
  })
  
  clusterExport(
    cl,
    varlist = c("evaluate_candidate_cluster", "MEE", "est_corr", "est_vgm", "vgm_plot", "Expvgm", "T.cnt")
  )
  
  # ------------------------------------------------------------------
  # Step 3: Fit Null Model & Estimate Spatiotemporal Correlation
  # ------------------------------------------------------------------
  message("Step 3: Fitting Null Model and estimating spatiotemporal correlation structure...")
  all.n.plus <- sum(response_mat > 0)
  P <- T.cnt
  q <- 1
  
  m0.0 <- MEE(
    response.mat  = response_mat,  
    covariate.mat = covariate_mat,  
    coordinate    = coordinate,  
    T.cnt         = T.cnt,  
    indep         = TRUE
  )
  
  m0.corr <- est_corr(
    epsilon.mat  = m0.0$epsilon,  
    response.mat = response_mat,  
    coordinate   = coordinate,  
    cutoff       = cutoff
  )
  
  sigma2.0       <- m0.corr$sigma2
  corr.set0      <- m0.corr$corr.set
  corr.plus.inv0 <- m0.corr$corr.ST.mat.plus.inv
  
  m0 <- MEE(
    response.mat  = response_mat,  
    covariate.mat = covariate_mat,  
    coordinate    = coordinate,  
    T.cnt         = T.cnt,  
    set.corr      = corr.set0,  
    sigma2        = sigma2.0,
    corr.plus.inv = corr.plus.inv0
  )
  
  SSE0.step1 <- m0$SSE
  
  # ------------------------------------------------------------------
  # Step 4: Step 1 Parallelized Candidate Scanning
  # ------------------------------------------------------------------
  message("Step 4: Executing parallelized candidate cluster scan (Step 1)...")
  sI <- splitIndices(nclust, length(cl))
  
  clusterExport(
    cl,
    varlist = c("response_mat", "covariate_mat", "coordinate", "Lambda", "cutoff", "corr.plus.inv0", "corr.set0", "T.cnt")
  )
  
  tmp1 <- clusterApply(
    cl,  
    x = sI,  
    fun = lapply,
    function(i) {
      evaluate_candidate_cluster(
        i             = i,  
        response.mat  = response_mat,  
        covariate.mat = covariate_mat,  
        coordinate    = coordinate,  
        Lambda        = Lambda,  
        cutoff        = cutoff,  
        corr.plus.inv = corr.plus.inv0,
        corr.set      = corr.set0,  
        T.cnt         = T.cnt
      )
    }
  )
  
  tmp.step1 <- do.call("c", tmp1)
  
  mtype.step1 <- sapply(tmp.step1, function(x) x$model)
  SSE1.step1  <- sapply(tmp.step1, function(x) x$SSE)
  SSE1.step1[SSE1.step1 > SSE0.step1] <- SSE0.step1  
  
  df1 <- P
  df2 <- all.n.plus - (q + 1) * Time - P
  
  F.step1 <- ((SSE0.step1 - SSE1.step1) / df1) / (SSE1.step1 / df2)
  p.step1 <- pf(q = F.step1, df1 = df1, df2 = df2, lower.tail = FALSE)
  p.step1[is.na(p.step1)] <- 1
  p.step1[mtype.step1 != "corr"] <- 1
  
  K0.step1 <- sum(p.step1 <= (0.05 / nclust), na.rm = TRUE)
  
  C.step1 <- Lambda[, order(p.step1)[1:K0.step1], drop = FALSE]
  K0.step1 <- ncol(C.step1)
  message(paste("-> Significant clusters identified in Step 1:", K0.step1))
  
  # ------------------------------------------------------------------
  # Step 5: Step 2 Cluster Regrouping
  # ------------------------------------------------------------------
  message("Step 5: Executing cluster regrouping (Step 2)...")
  C.step2 <- C.step1[, 1, drop = FALSE]
  
  if (K0.step1 > 1) {
    for (i in 2:K0.step1) {
      new.clust <- C.step1[, i, drop = FALSE]
      
      if (all((1 - (rowSums(C.step2) > 0)) * new.clust == 0)) next
      if (all((rowSums(C.step2) * new.clust) == 0)) {
        C.step2 <- cbind(C.step2, new.clust)
        next
      }
      
      sel.loc <- max(which(colSums(C.step2 * c(new.clust)) > 0))
      one.clust <- C.step2[, sel.loc, drop = FALSE]
      new.clust <- (1 - (rowSums(C.step2) > 0)) * new.clust
      
      if (sum(new.clust) <= 3) {
        C.step2[, sel.loc] <- C.step2[, sel.loc] + new.clust
        next
      }
      
      two.clust <- cbind(one.clust, new.clust)
      
      m1.0 <- MEE(response.mat = response_mat, covariate.mat = covariate_mat, clust = one.clust, coordinate = coordinate, T.cnt = T.cnt, indep = TRUE)
      sigma2.1 <- var(c(m1.0$epsilon), na.rm = TRUE)
      m1 <- MEE(response.mat = response_mat, covariate.mat = covariate_mat, clust = one.clust, coordinate = coordinate, T.cnt = T.cnt, set.corr = corr.set0, sigma2 = sigma2.1, corr.plus.inv = corr.plus.inv0)
      
      m2.0 <- MEE(response.mat = response_mat, covariate.mat = covariate_mat, clust = two.clust, coordinate = coordinate, T.cnt = T.cnt, indep = TRUE)
      sigma2.2 <- var(c(m2.0$epsilon), na.rm = TRUE)
      m2 <- MEE(response.mat = response_mat, covariate.mat = covariate_mat, clust = two.clust, coordinate = coordinate, T.cnt = T.cnt, set.corr = corr.set0, sigma2 = sigma2.2, corr.plus.inv = corr.plus.inv0)
      
      SSE1.step2 <- m1$SSE
      SSE2.step2 <- min(m2$SSE, SSE1.step2)
      
      df1.step2 <- P
      df2.step2 <- all.n.plus - (q + 1) * Time - 2 * P
      
      F.step2 <- ((SSE1.step2 - SSE2.step2) / df1.step2) / (SSE2.step2 / df2.step2)
      p.step2 <- pf(q = F.step2, df1 = df1.step2, df2 = df2.step2, lower.tail = FALSE)
      
      if (p.step2 < 0.05) {
        C.step2 <- cbind(C.step2, new.clust)
      } else {
        C.step2[, sel.loc] <- C.step2[, sel.loc] + new.clust
      }
    }
  }
  
  # ------------------------------------------------------------------
  # Step 6: Final Model Refinement, Shapefile Update & Visualization
  # ------------------------------------------------------------------
  message("Step 6: Fitting final model, updating shapefile, and generating visualization...")
  m3.0 <- MEE(response.mat = response_mat, covariate.mat = covariate_mat, clust = C.step2, coordinate = coordinate, T.cnt = T.cnt, indep = TRUE)
  m3.corr <- est_corr(epsilon.mat = m3.0$epsilon, response.mat = response_mat, coordinate = coordinate, cutoff = cutoff)
  
  sigma2.3       <- m3.corr$sigma2
  corr.set3      <- m3.corr$corr.set
  corr.plus.inv3 <- m3.corr$corr.ST.mat.plus.inv
  
  final.model <- MEE(
    response.mat  = response_mat,  
    covariate.mat = covariate_mat,  
    coordinate    = coordinate,  
    T.cnt         = T.cnt,
    set.corr      = corr.set3, 
    sigma2        = sigma2.3,
    corr.plus.inv = corr.plus.inv3
  )
  
  stopCluster(cl)
  
  # Attach cluster results to the shapefile object
  K_final <- ncol(C.step2)
  cluster_code <- c(C.step2 %*% c(1:K_final))
  updated_shp <- shp.file
  updated_shp$cluster_id    <- cluster_code
  updated_shp$cluster_group <- ifelse(cluster_code > 0, paste0("C", cluster_code), "Non-Cluster")
  
  # Generate Map Plot
  fig.final <- clusters_plot(clust = C.step2, shp.file = updated_shp)
  fig.final <- fig.final + theme(text = element_text(family = "serif")) +
    ggtitle(paste0("Final Clusters: ", K_final))
  
  t1 <- Sys.time()
  message(paste("Total execution time:", round(difftime(t1, t0, units = "mins"), 2), "minutes"))
  
  list(
    final.model    = final.model,
    final.clusters = C.step2,
    updated.shp    = updated_shp,
    plot           = fig.final
  )
}