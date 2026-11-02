#' Mixture Estimation Equations (MEE) Model for Spatiotemporal Data
#'
#' @description Fits a spatiotemporal mixture regression model for zero-inflated or mixed 
#' responses (comprising zeros and continuous positive values) and evaluates cluster significance 
#' by combining spatiotemporal correlation structures with cluster indicator variables via 
#' iterative reweighted least squares (IRLS) and score derivatives.
#'
#' @param response.mat Matrix of response variables (rows: spatial locations, columns: time points).
#' @param covariate.mat Matrix or array of covariates (comprising zeros and continuous positive values).
#' @param clust Matrix, candidate cluster indicator matrix (defaults to \code{NULL}).
#' @param coordinate Spatial coordinate matrix.
#' @param tol Convergence threshold (defaults to \code{10^(-4)}).
#' @param max.iter Maximum number of iterations (defaults to 20).
#' @param set.corr Vector of pre-specified spatiotemporal correlation parameters.
#' @param corr.plus.inv Inverse spatiotemporal correlation matrix for positive observations.
#' @param indep Logical value indicating whether spatial and temporal independence is assumed (defaults to \code{FALSE}).
#' @param sigma2 Initial error variance value (defaults to \code{NULL}).
#' @param T.cnt Temporal trend degree (1 for intercept, 2 for linear, 3 for quadratic; defaults to 3).
#'
#' @return A list containing convergence status, iteration count, coefficient tables, standard errors, 
#' sandwich covariance matrices, fitted values, residuals, SSE, residual variance, and temporal trend components.
#' 
#' @export
MEE <- function(response.mat, covariate.mat, clust = NULL, coordinate, 
                tol = 10^(-4), max.iter = 20, set.corr = NULL, 
                corr.plus.inv = NULL, indep = FALSE, sigma2 = NULL, T.cnt = 3) {
  
  Time <- ncol(response.mat)
  n <- nrow(response.mat)
  q <- 1
  na.loc <- NULL
  L2 <- as.matrix(dist(coordinate))
  
  X.list <- lapply(1:Time, function(t) cbind(1, covariate.mat[, t]))
  X_block <- bdiag(X.list)
  
  if(!is.null(clust)) {
    K <- ncol(clust)
    big.C <- kronecker(diag(Time), clust)
    X <- cbind(X_block, big.C)
    R_C <- clust %*% t(clust)
  } else {
    X <- X_block
    R_C <- matrix(0, n, n)
  }
  
  Y <- c(response.mat)
  part.plus <- Y > 0
  n.plus <- sum(part.plus)
  Y.plus <- Y[part.plus]
  X.plus <- X[part.plus, , drop = FALSE]
  Big.Psi <- 1 - part.plus
  
  beta_len <- Time * (q + 1)
  X_0 <- as.matrix(X)
  X_0[ , -c(1:beta_len)] <- 0
  X_0 <- Matrix(X_0, sparse = TRUE)
  
  # Obtain initial regression coefficients using a traditional Gaussian linear model on positive subsets
  fit0 <- glm(Y.plus ~ 0 + as.matrix(X.plus), family = gaussian)
  old.Beta <- fit0$coefficients
  old.Beta[is.na(old.Beta)] <- 0
  
  if(is.null(sigma2)) {
    sigma2 <- var(fit0$residuals)
  }
  sigma.hat <- sqrt(sigma2)
  
  if(indep){
    R_T <- diag(Time)
    R_S <- diag(n)
    corr.coef <- c(c0 = NA, c1 = NA, range = NA, ar = NA)
  } else {
    corr.coef <- set.corr
    if(is.null(corr.plus.inv)){
      stop("Error: 'corr.plus.inv' must be provided when indep = FALSE.")
    }
    covY.plus.inv <- (sigma.hat ^ (-2)) * corr.plus.inv
  }
  
  iter <- 0
  
  # Enter iterative optimization loop (Newton-Raphson / Fisher Scoring algorithm)
  repeat {
    iter <- iter + 1
    
    Y.hat.plus <- as.numeric(X.plus %*% old.Beta)
    epsilon.plus <- Y.plus - Y.hat.plus
    
    Y.hat <- as.numeric(X %*% old.Beta)
    epsilon <- Y - Y.hat
    
    all.epsilon <- Y - as.numeric(X %*% old.Beta)
    all.epsilon[Y == 0] <- 0
    all.epsilon.mat <- matrix(all.epsilon, nrow = n)
    epsilon.T <- apply(all.epsilon.mat, 2, mean)
    epsilon.S <- apply(all.epsilon.mat, 1, mean)
    
    if(indep){
      sigma.hat <- sd(c(all.epsilon))
      covY.plus.inv <- (sigma.hat ^ (-2)) * Diagonal(n = sum(part.plus))
    }
    
    Z.X <- covY.plus.inv %*% X.plus
    
    phi_0 <- pi / (sqrt(3) * sigma.hat)
    theta.hat <- as.numeric((1 + exp(phi_0 * X %*% old.Beta))^(-1)) 
    D.theta <- Diagonal(x = theta.hat * (1 - theta.hat))
    
    # Compute Score function (L) and approximate Hessian matrix (L.dot)
    L <- phi_0 * t(X_0) %*% (theta.hat - Big.Psi) + 
      t(X.plus) %*% covY.plus.inv %*% (epsilon.plus)
    L.dot <- - ((phi_0 ^ 2) * t(X_0) %*% D.theta %*% X_0 + 
                  t(X.plus) %*% Z.X)
    
    new.Beta <- try(old.Beta - solve(L.dot, L), silent = TRUE)
    
    if(inherits(new.Beta, "try-error")) {
      new.Beta <- old.Beta - ginv(as.matrix(L.dot)) %*% as.matrix(L)
    }
    
    if(sqrt(sum((new.Beta - old.Beta)^2)) < tol) break 
    if(iter > max.iter) {
      warning("Warning: MEE model failed to converge within maximum iterations.")
      break
    }
    old.Beta <- new.Beta
  }
  
  # Calculate robust sandwich covariance matrix and standard errors
  Sigma.mat1 <- t(X.plus) %*% Z.X
  Sigma.mat0 <- (sigma.hat ^ (-2)) * t(X) %*% X
  sand.Sigma_Matrix <- try(solve(Sigma.mat0) %*% Sigma.mat1 %*% solve(Sigma.mat0), silent = TRUE)
  
  if(inherits(sand.Sigma_Matrix, "try-error")) {
    sand.Sigma_Matrix <- try(ginv(as.matrix(Sigma.mat0)) %*% Sigma.mat1 %*% ginv(as.matrix(Sigma.mat0)), silent = TRUE)
  }
  
  if(inherits(sand.Sigma_Matrix, "try-error")) {
    eps <- 1e-6
    Sigma.reg <- Sigma.mat0 + diag(eps, nrow(Sigma.mat0))
    sand.Sigma_Matrix <- solve(Sigma.reg) %*% Sigma.mat1 %*% solve(Sigma.reg)
  }
  
  Sigma_Matrix <- try(solve(Sigma.mat1), silent = TRUE)
  if(inherits(Sigma_Matrix, "try-error")) {
    Sigma_Matrix <- try(ginv(as.matrix(Sigma.mat1)), silent = TRUE)
  }
  if(inherits(Sigma_Matrix, "try-error")) {
    eps <- 1e-6
    Sigma.reg <- Sigma.mat1 + diag(eps, nrow(Sigma.reg))
    Sigma_Matrix <- solve(Sigma.reg)
  }
  
  res.coef <- as.numeric(new.Beta)
  res.sd <- sqrt(diag(Sigma_Matrix))
  res.zval <- res.coef / res.sd
  res.pval <- 2 * (1 - pnorm(abs(res.zval)))
  
  if (!is.null(clust)) {
    beta_len <- Time * (q + 1)
    
    xi_coefs <- res.coef[(beta_len + 1):length(res.coef)]
    beta_coef <- res.coef[1:beta_len]
    xi_mat <- Phi <- matrix(xi_coefs, nrow = Time, ncol = K, byrow = TRUE)
    
    # Build temporal trend matrix D
    if (T.cnt == 1) {
      D <- matrix(1, nrow = Time, ncol = 1)
    } else if (T.cnt == 2) {
      D <- cbind(1, 1:Time)
    } else {
      D <- cbind(1, 1:Time, (1:Time)^2)
    }
    
    # OLS and GLS temporal trend fitting
    phi.ols <- solve(t(D) %*% D) %*% t(D) %*% Phi
    eta.hat <- Phi - D %*% phi.ols
    
    acf.list <- vector("list", K)
    for (k in 1:K) {
      acf.list[[k]] <- acf(eta.hat[, k], lag.max = Time - 1, plot = FALSE)
    }
    
    Sigma.eta <- vector("list", K)
    sigma2.eta <- numeric(K)
    rho.eta <- matrix(NA, nrow = K, ncol = Time)
    phi.gls <- matrix(NA, nrow = ncol(D), ncol = K)
    
    for (k in 1:K) {
      rho.eta.k <- as.numeric(acf.list[[k]]$acf)
      rho.eta[k, ] <- rho.eta.k
      sigma2.eta[k] <- var(eta.hat[, k])
      gamma.eta <- sigma2.eta[k] * rho.eta.k
      Sigma.eta[[k]] <- toeplitz(gamma.eta)
      
      Sigma.eta.inv <- solve(Sigma.eta[[k]])
      phi.gls[, k] <- solve(t(D) %*% Sigma.eta.inv %*% D) %*% t(D) %*% Sigma.eta.inv %*% Phi[, k]
    }
    
    Phi.gls.fitted <- D %*% phi.gls
    Phi.gls.fitted.vec <- as.vector(t(Phi.gls.fitted))
    res.fitted <- as.numeric(X %*% c(beta_coef, Phi.gls.fitted.vec))
    res.T.model <- Phi.gls.fitted.vec
    phi.mat <- phi.gls
    
  } else {
    res.fitted <- as.numeric(X %*% res.coef)
    phi.mat <- D <- NULL
    sigma2.eta <- rho.eta <- NULL
  }
  
  res.epsilon <- res.epsilon.star <- Y - res.fitted
  res.epsilon.plus <- res.epsilon[part.plus]
  res.epsilon[!part.plus] <- 0
  res.epsilon <- matrix(res.epsilon, nrow = n, ncol = Time)
  res.epsilon.star <- matrix(res.epsilon.star, nrow = n, ncol = Time)
  
  res.SSE <- sum((res.epsilon.plus)^2)
  
  res.table <- res.table2 <- data.frame(
    coefficients = res.coef,
    standard.errors = res.sd,
    z.values = res.zval,
    p.values = res.pval
  )
  
  if(!is.null(names(corr.coef)) && length(corr.coef) == 4) {
    names(corr.coef) <- c("c0", "c1", "range", "ar")
  }
  
  list(
    converge          = iter <= max.iter, 
    iter              = iter, 
    table             = res.table, 
    table2            = res.table2, 
    sand.Sigma_Matrix = sand.Sigma_Matrix, 
    Sigma_Matrix      = Sigma_Matrix,
    fitted.values     = res.fitted, 
    epsilon           = res.epsilon, 
    epsilon.star      = res.epsilon.star,
    SSE               = res.SSE, 
    sigma2            = sigma.hat ^ 2, 
    corr.coef         = corr.coef,
    sigma2.eta        = sigma2.eta, 
    rho.eta           = rho.eta, 
    phi.mat           = phi.mat, 
    D                 = D
  )
}