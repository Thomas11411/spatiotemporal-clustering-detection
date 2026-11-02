#' Estimate Spatial Semivariance and Correlation Matrix
#'
#' @description Computes empirical semivariance from residual vectors and spatial coordinates, 
#' fits an exponential variogram model, and derives the spatial correlation matrix.
#'
#' @param epsilon Numeric vector of spatial residuals.
#' @param coordinate Matrix or data frame of 2D spatial coordinates (x, y).
#' @param cutoff Numeric value, the maximum distance cutoff for semivariance calculation (defaults to \code{NULL}).
#'
#' @return A list containing model existence status (\code{exist}), spatial correlation matrix (\code{corr}), 
#' nugget (\code{c0}), partial sill (\code{c1}), range (\code{range}), plot object (\code{p}), and fitted models.
#' 
#' @export
est_vgm <- function(epsilon, coordinate, cutoff = NULL) {
  
  data <- data.frame(e = epsilon, x = coordinate[, 1], y = coordinate[, 2])
  
  # Calculate empirical semivariogram
  if(is.null(cutoff)) {
    gamma.vgm <- variogram(e ~ 1, locations = ~ x + y, data = data)
  } else {
    gamma.vgm <- variogram(e ~ 1, locations = ~ x + y, 
                           cutoff = cutoff, data = data)
  }
  
  # Attempt non-linear least squares fitting using an Exponential variogram model
  fit.vgm <-
    tryCatch( 
      {fit.variogram(gamma.vgm, model = vgm("Exp"), fit.method = 2)},
      warning = function(w) { rho.key = 1 },
      error = function(e) { rho.key = 2 }
    )
  
  if(inherits(fit.vgm, "variogramModel")){
    
    # If the fitted range is less than or equal to 0, fitting fails and defaults to a diagonal identity matrix (independence)
    if(fit.vgm[2, "range"] <= 0){
      exist <- FALSE
      corr <- diag(1, length(epsilon))
      c0 <- NA
      c1 <- NA
      range <- NA
      sigma2 <- var(epsilon)
      p <- vgm_plot(v = gamma.vgm)
    } else {
      exist <- TRUE
      L2 <- as.matrix(dist(coordinate))
      c0 <- fit.vgm[1, "psill"]
      c1 <- fit.vgm[2, "psill"]
      range <- fit.vgm[2, "range"]
      
      # Compute spatial semivariance matrix based on the exponential variogram model
      gammahat.vgm <- apply(L2, c(1, 2), Expvgm, c0 = c0, c1 = c1, r = range)
      # Convert semivariance to spatial correlation matrix: Corr = 1 - (Gamma / Sill)
      corr <- 1 - gammahat.vgm / (c0 + c1)
      p <- vgm_plot(v = gamma.vgm, m = fit.vgm)
    }
    
  } else {
    exist <- FALSE
    corr <- diag(1, length(epsilon))
    c0 <- NA
    c1 <- NA
    range <- NA
    sigma2 <- var(epsilon)
    p <- vgm_plot(v = gamma.vgm)
  }
  
  list(exist = exist, corr = corr, c0 = c0, c1 = c1, range = range, p = p, gamma.vgm = gamma.vgm, fit.vgm = fit.vgm)
}
