#' Two-Stage Model Fitting and SSE Evaluation for a Single Candidate Cluster
#'
#' @description Fits a regression model for a specified candidate cluster in two stages to compute its Residual Sum of Squares (SSE):
#' Stage 1 fits an independence MEE model to obtain stable initial values, baseline variance, and independence SSE (\code{SSE.indep}).
#' Stage 2 incorporates the spatiotemporal correlation structure to compute the final model SSE (\code{SSE}), falling back safely 
#' to the independence model if convergence fails.
#'
#' @param i Numeric index of the candidate cluster.
#' @param response.mat Response matrix.
#' @param covariate.mat Covariate matrix.
#' @param Lambda Indicator matrix for all candidate clusters.
#' @param coordinate Spatial coordinate matrix.
#' @param corr.plus.inv Inverse spatiotemporal correlation matrix.
#' @param corr.set Correlation parameters.
#' @param cutoff Spatial distance cutoff.
#' @param T.cnt Temporal trend degree.
#'
#' @return A list containing the final model type (\code{corr} or \code{indep}), sandwich covariance matrix, 
#' regression coefficients, variance \code{sigma2}, correlation parameters, final SSE (\code{SSE}), 
#' independence SSE (\code{SSE.indep}), and execution time in seconds.
#' 
#' @export
evaluate_candidate_cluster <- function(i, response.mat, covariate.mat, Lambda, coordinate, 
                                       corr.plus.inv, corr.set, cutoff, T.cnt) {
  
  ti.0 <- Sys.time()
  Lambda_i <- Lambda[, i, drop = FALSE]
  
  # ------------------------------------------------------------------
  # Stage 1: Independence Model Fitting (indep = TRUE) for initial values and baseline SSE
  # ------------------------------------------------------------------
  m1.0 <- try(
    MEE(
      response.mat  = response.mat,     
      covariate.mat = covariate.mat,     
      clust         = Lambda_i,         
      coordinate    = coordinate, 
      T.cnt         = T.cnt, 
      indep         = TRUE
    ), 
    silent = TRUE
  )
  
  # If Stage 1 fails or fails to converge, return default null/NA values immediately
  if (inherits(m1.0, "try-error") || !isTRUE(m1.0$converge)) {
    return(list(
      model      = NA,
      sand.S.M  = NULL,
      S.M       = NULL,
      coef      = NULL,
      sigma2    = NA,
      corr.set  = rep(NA, 4),
      SSE       = NA,
      SSE.indep = NA,
      time_sec  = as.numeric(Sys.time() - ti.0, units = "secs")
    ))
  }
  
  SSE.indep <- m1.0$SSE
  sigma2.1  <- var(c(m1.0$epsilon), na.rm = TRUE)
  
  # ------------------------------------------------------------------
  # Stage 2: MEE Model Fitting Incorporating Spatiotemporal Correlation Structure
  # ------------------------------------------------------------------
  m1 <- try(
    MEE(
      response.mat  = response.mat,     
      covariate.mat = covariate.mat,     
      clust         = Lambda_i,         
      coordinate    = coordinate, 
      T.cnt         = T.cnt,
      set.corr      = corr.set,     # Pass correct spatiotemporal correlation parameters
      sigma2        = sigma2.1,
      corr.plus.inv = corr.plus.inv
    ), 
    silent = TRUE
  )
  
  # ------------------------------------------------------------------
  # Check Stage 2 Convergence: Use correlated model SSE if successful, fallback to independence SSE otherwise
  # ------------------------------------------------------------------
  if (!inherits(m1, "try-error") && isTRUE(m1$converge)) {
    model     <- "corr"
    sand.S.M <- m1$sand.Sigma_Matrix
    S.M      <- m1$Sigma_Matrix
    coef     <- m1$table[, 1]
    sigma2   <- sigma2.1
    res.corr <- corr.set
    SSE      <- m1$SSE
  } else {
    model     <- "indep"
    sand.S.M <- m1.0$sand.Sigma_Matrix
    S.M      <- m1.0$Sigma_Matrix
    coef     <- m1.0$table[, 1]
    sigma2   <- m1.0$sigma2
    res.corr <- rep(NA, 4)
    SSE      <- SSE.indep
  }
  
  ti.1 <- Sys.time()
  
  list(
    model       = model,
    sand.S.M  = sand.S.M, 
    S.M       = S.M,
    coef      = coef, 
    sigma2    = sigma2, 
    corr.set  = res.corr,
    SSE       = SSE, 
    SSE.indep = SSE.indep,
    time_sec  = as.numeric(ti.1 - ti.0, units = "secs")
  )
}