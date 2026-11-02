#' Estimate Spatiotemporal Joint Correlation Structure (Kronecker Product Structure)
#'
#' @description Independently estimates spatial correlation (via variogram models) and temporal correlation 
#' (via AR(1) models) from residual matrices, constructs the joint spatiotemporal correlation matrix using a 
#' Kronecker product, and computes its inverse securely for Generalized Estimating Equations (GEE).
#'
#' @param epsilon.mat Matrix of spatiotemporal residuals (rows represent spatial locations, columns represent time points).
#' @param response.mat Matrix of response variables.
#' @param coordinate Spatial coordinate matrix.
#' @param cutoff Numeric cutoff distance for spatial variogram.
#'
#' @return A list containing total variance \code{sigma2}, correlation parameter set \code{corr.set}, 
#' spatial and temporal correlation matrices, and the subset-filtered inverse spatiotemporal matrix \code{corr.ST.mat.plus.inv}.
#' 
#' @export
est_corr <- function(epsilon.mat, response.mat, coordinate, cutoff) {
  
  # Identify positive response observations for zero-inflation or subset filtering
  part.plus <- c(response.mat > 0)
  
  epsilon.S <- apply(epsilon.mat, 1, mean)
  epsilon.T <- apply(epsilon.mat, 2, mean)
  
  n <- length(epsilon.S)
  Time <- ncol(epsilon.mat)
  
  # Estimate total spatiotemporal variance
  sigma2 <- var(c(epsilon.mat), na.rm = TRUE)
  
  # ------------------------------------------------------------
  # 1. Estimate Spatial Correlation
  # ------------------------------------------------------------
  vgm.fit <- est_vgm(
    epsilon = epsilon.S,
    coordinate = coordinate,
    cutoff = cutoff
  )
  
  vgm.coef <- c(
    vgm.fit$c0,
    vgm.fit$c1,
    vgm.fit$range
  )
  
  corr.S.mat <- vgm.fit$corr
  
  # ------------------------------------------------------------
  # 2. Estimate Temporal Correlation (AR(1) Model)
  # ------------------------------------------------------------
  corr.T <- ar(
    epsilon.T,
    order.max = 1,
    aic = FALSE
  )$ar
  
  # Ensure the AR(1) coefficient stays within the stationarity open interval (0, 1)
  corr.T <- ifelse(corr.T > 0 && corr.T < 1, corr.T, 0)
  time_diff <- abs(outer(c(1:Time), c(1:Time), "-"))
  corr.T.mat <- corr.T ^ time_diff
  
  # ------------------------------------------------------------
  # 3. Aggregate Correlation Parameters and Compute Kronecker Product
  # ------------------------------------------------------------
  corr.set <- c(vgm.coef, corr.T)
  names(corr.set) <- c("c0", "c1", "range", "ar")
  
  # Build the full spatiotemporal correlation matrix via Kronecker product
  corr.ST.mat <- kronecker(corr.T.mat, corr.S.mat)
  corr.ST.mat.plus <- corr.ST.mat[part.plus, part.plus, drop = FALSE]
  
  # Safe matrix inversion handling: fallback to generalized pseudo-inverse (ginv) if nearly singular
  corr.ST.mat.plus.inv <- tryCatch(
    solve(corr.ST.mat.plus),
    error = function(e) {
      warning("Matrix is singular; falling back to generalized inverse (ginv)")
      MASS::ginv(corr.ST.mat.plus)
    }
  )
  
  list(
    sigma2 = sigma2,
    corr.set = corr.set,
    corr.S.mat = corr.S.mat,
    corr.T.mat = corr.T.mat,
    corr.ST.mat.plus.inv = corr.ST.mat.plus.inv,
    p = vgm.fit$p
  )
}