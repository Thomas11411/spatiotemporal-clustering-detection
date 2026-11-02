#' Exponential Variogram Model
#'
#' @description Computes the theoretical semivariance at a given spatial distance \code{h}.
#'
#' @param h Numeric value or vector representing spatial distance.
#' @param c0 Numeric value, the nugget effect.
#' @param c1 Numeric value, the partial sill.
#' @param r Numeric value, the range parameter of the variogram.
#'
#' @return Returns the theoretical semivariance value \code{gamma} corresponding to distance \code{h}.
#' 
#' @export
Expvgm <- function(h, c0, c1, r) {
  if ( h == 0 ) { 
    gamma <- 0 
  } else { 
    gamma <- c0 + c1 * (1 - exp(-h / r)) 
  }
  return(gamma) 
}