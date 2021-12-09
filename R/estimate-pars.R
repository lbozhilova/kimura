#' Estimate Kimura distribution parameters
#'
#' Given a vector of heteroplasmies, estimated the Kimura distribution
#' parameters.
#'
#' @param h numeric: heteroplasmy fractions
#'
#' @return A named vector of length two.
#' @export
#'
#' @examples
#' h = stats::runif(10)
#' estimate_parameters(h)
estimate_parameters <- function(h) {
  p <- mean(h)
  b <- 1 - stats::var(h) / (p * (1 - p))
  c("p" = p, "b" = b)
}

#' The Kimura distribution
#'
#' Get complete CDF for Kimura(p, b).
#'
#' @inheritParams .phi
#'
#' @return list of length 2, x and y
#' @export
#'
#' @examples
#' pkimura_full(0.3, 0.5)
pkimura_full <- function(p, b) {
  x <- seq(0, 1, 1e-4)
  pdf_x <- sapply(x, function(x) dkimura(x, p, b))
  cdf_x <- cumsum(pdf_x * c(1, rep(1e-4, 10000)))
  cdf_x[10001] <- 1
  list("x" = x ,
       "y" = cdf_x)
}
