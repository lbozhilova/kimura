#' Monte Carlo Kolmogorov-Smirnov test for genetic drift
#'
#' Given  heteroplasmy fractions \code{h}, test whether their distribution is
#' consistent with random genetic drift.
#'
#' @inheritParams estimate_parameters
#' @param num_MC number of Monte Carlo simulations.
#' @param round logical, indicates whether heteroplasmy fractions are rounded to
#'   two significant digits.
#'
#' @return object of class htest
#' @export
#'
#' @examples
#' h <- seq(0.2, 0.5, 0.01)
#' test_kimura(h)
test_kimura <- function(h, num_MC = 1000, round = TRUE) {
  pars <- estimate_parameters(h)
  p <- unname(pars[1])
  b <- unname(pars[2])
  cdf_kimura <- .pkimura_full(p, b)
  h_mat <- cbind(h, matrix(rkimura(num_MC * length(h), p, b), ncol = num_MC))
  if (round)
    h_mat <- round(h_mat, 2)
  get_ks_statistic <- function(h) {
    ecdf_h <- stats::ecdf(h)(seq(0, 1, 1e-4))
    max(abs(ecdf_h - cdf_kimura))
  }
  D <- unlist(apply(h_mat, 2, get_ks_statistic))
  D <- unname(D)
  output <- list(
    statistic = c("D" = D[1], "p" = p, "b" = b),
    p.value = unname(sum(D >= D[1]) / length(D)),
    alternative = "one-sided",
    method = "Monte Carlo Kolmogorov-Smirnov",
    data.name = paste0(deparse(substitute(h)),
                       " and Kimura(", round(p, 4), ", ", round(b, 4), ")"),
    D_sample = D)
  class(output) <- "htest"
  output
}
