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

