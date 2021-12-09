#' Fit a Kimura distirbution to heteroplasmy data
#'
#' @inheritParams  estimate_parameters
#'
#' @return List of length 3 = p, b, and the estomate
#' @export
#'
#' @examples
#' h <- seq(0.2, 0.5, 0.01)
#' fit_kimura(h)
fit_kimura <- function(h) {
  pars <- estimate_parameters(h)
  p <- pars[1]
  b <- pars[2]
  cdf_lt <- pkimura_full(p, b)
  list(
    "p" = p,
    "b" = b,
    "cdf" = cdf_lt)
}

#' Get the KS test statistic between a heteroplasmy CDF and a Kimura CDF
#'
#' @inheritParams  estimate_parameters
#' @param cdf_lt list of length 2, see pkimura_full
#'
#' @return Positive D
#' @export
#'
#' @examples
#' h <- seq(0.2, 0.5, 0.01)
#' cdf_lt <- fit_kimura(h)$cdf
#' .get_ks_stat(h, cdf_lt)
.get_ks_stat <- function(h, cdf_lt) {
  ecdf_h <- stats::ecdf(h)
  cdf_h_eval <- ecdf_h(cdf_lt$x)
  max(abs(cdf_h_eval - cdf_lt$y))
}

#' Sample heteroplasmies from a given Kimura CDF
#'
#' @inheritParams .get_ks_stat
#' @inheritParams rkimura
#' @param round logical; whether to round to 2 decimal places
#'
#' @return vector length n
#' @export
#'
#' @examples
#' h <- seq(0.2, 0.5, 0.01)
#' cdf_lt <- fit_kimura(h)
#' .sample_h(length(h), cdf_lt)
.sample_h <- function(n, cdf_lt, round = TRUE) {
  p <- cdf_lt$p
  b <- cdf_lt$b
  spl <- sample(c(0, 0.5, 1), n, replace = TRUE,
                prob = c(.f0(p, b), 1 - .f0(p, b) - .f1(p, b), .f1(p, b)))
  idx <- which(spl == 0.5)
  m <- length(idx)
  spl[idx] <- stats::approx(cdf_lt$cdf$y,
                            cdf_lt$cdf$x,
                            stats::runif(m, min = cdf_lt$cdf$y[1]))$y
  if (round)
    spl <- round(spl, 2)
  spl
}

#' Monte Carlo KS test for genetic drift
#'
#' @param h heteroplasmy
#' @param num_MC number MC simulations
#' @param round whether the het is rounded
#'
#' @return object of class htest
#' @export
#'
#' @examples
#' h <- seq(0.2, 0.5, 0.01)
#' test_kimura(h)
test_kimura <- function(h, num_MC = 1000, round = TRUE) {
  fitted <- fit_kimura(h)
  s_h <- matrix(
    .sample_h(num_MC * length(h), fitted, round),
    ncol = num_MC)
  ks_obs <- .get_ks_stat(h, fitted$cdf)
  ks_spl <- apply(s_h, 2, function(sh) .get_ks_stat(sh, fitted$cdf))
  ks_spl <- unlist(ks_spl)
  output <- list(
    statistic = c("D" = ks_obs, "p" = fitted$p, "b" = fitted$b),
    p.value = (1 + sum(ks_spl >= ks_obs)) / (num_MC + 1),
    alternative = "one-sided",
    method = "Monte Carlo Kolmogorov-Smirnov",
    data.name = paste0("h and Kimura", fitted$p, ", ", fitted$b, ")"),
    D_sample = ks_spl)
  class(output) <- "htest"
  output
}
