#' Estimate Kimura distribution parameters
#'
#' Given a vector of heteroplasmy fractions, estimate the Kimura distribution
#' parameters.
#'
#' @param h numeric vector in \eqn{[0; 1]} of heteroplasmy fractions.
#'
#' @return A named vector of length two with elements \code{p} and \code{b}.
#'   These are the estimated parameters of a \eqn{Kimura(p, b)} distribution to
#'   fit to the data.
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

#' The full Kimura distribution CDF
#'
#' Get the CDF for \code{Kimura(p, b)} evaluated at regular intervals of length
#' 1e-04. Note: this is an auxiliary function.
#'
#' @inheritParams .phi
#'
#' @return Numeric vector of length 10001, corresponding to the CDF for
#'   \code{Kimura(p, b)}, approximated at \code{x = seq(0, 1, 1e-4)}.
#'
#' @references Wonnapinij, Passorn, Patrick F. Chinnery, and David C. Samuels.
#'   "The distribution of mitochondrial DNA heteroplasmy due to random genetic
#'   drift." The American Journal of Human Genetics 83.5 (2008): 582-593.
#' @export
#'
#' @examples
#' .pkimura_full(0.3, 0.5)
.pkimura_full <- function(p, b) {
  x <- seq(0, 1, 1e-4)
  pdf_x <- dkimura(x, p, b)
  cdf_x <- numeric(length(x))
  cdf_x[1] <- pdf_x[1]
  cdf_x[2] <- pdf_x[1] + 0.5 * 1e-4 * pdf_x[2]
  for (i in 3:(length(x)-1))
    cdf_x[i] <- cdf_x[i-1] + 0.5 * 1e-4 * (pdf_x[i-1] + pdf_x[i+1])
  cdf_x[length(x)] <- 1
  cdf_x
}
