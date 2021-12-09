#' Estimate F(i-1, i+2, 2, x)
#'
#' Estimates the hypergeometric function F(i-1, i+2, 2, x) needed for
#' calculating the Kimura distribution.
#'
#' @param i positive integer: the index in the infinite sum of the Kimura distribution.
#' @param x number in (0; 1): the point at which F() is calculated.
#'
#' @return The value of F(i-1, i+2, 2, x), estimated recursively.
#' @export
#'
#' @examples
#' .hypgeo(5, .3)
.hypgeo <- function(i, x) {
  a <- 1 - i
  b <- i + 2
  c <- 2
  f0 <- 1
  f1 <- 1 - (b * x) / c
  if (i == 1)
    return(f0)
  if (i == 2)
    return(f1)
  f <- 0
  for (j in 2:(i-1)) {
    aa <- 1 - j
    f <- (aa * (1 - x) * (f1 - f0) + (aa + b * x - c) * f1) / (aa - c)
    f0 <- f1
    f1 <- f
  }
  f
}

#' Probability of allele loss
#'
#' Estimates the probability that an allele is lost due to random genetic drift.
#'
#' @param p number in (0; 1): initial heteroplasmy
#' @param b number in (0; 1): time / rate parameter
#'
#' @return The probability the allele is lost due to genetic drift.
#' @export
#'
#' @examples
#' .f0(0.3, 0.5)
.f0 <- function(p, b) {
  q <- 1 - p
  N <- 250
  sum_terms <- numeric(N)
  i_coef <- p * q * (2 * (1:N) + 1) * c(-1, 1)
  b_coef <- b^choose(2:(N + 1), 2)
  q_hg <- .hypgeo(1, q)
  sum_terms[1] <- q_hg * i_coef[1] * b_coef[1]
  for (i in 2:N) {
    q_hg <- .hypgeo(i, q)
    sum_terms[i] <- q_hg * i_coef[i] * b_coef[i]
    if (abs(sum_terms[i]-sum_terms[i-1]) < 1e-5)
      break
  }
  if(abs(sum_terms[N - 1] - sum_terms[N]) > 1e-5)
    warning("Series not yet converged at N = 250.")
  q + sum(sum_terms[1:i])
}

#' Probability of fixing allele
#'
#' Estimates the probability that an allele is fixed (i.e. the wildtype is lost)
#' due to random genetic drift.
#'
#' @inheritParams .f0
#'
#' @return The probability the allele is fixed due to genetic drift.
#' @export
#'
#' @examples
#' .f1(0.3, 0.5)
.f1 <- function(p, b) {
  q <- 1 - p
  N <- 250
  sum_terms <- numeric(N)
  i_coef <- p * q * (2 * (1:N) + 1) * c(-1, 1)
  b_coef <- b^choose(2:(N + 1), 2)
  p_hg <- .hypgeo(1, p)
  sum_terms[1] <- p_hg * i_coef[1] * b_coef[1]
  for (i in 2:N) {
    p_hg <- .hypgeo(i, p)
    sum_terms[i] <- p_hg * i_coef[i] * b_coef[i]
    if (abs(sum_terms[i]-sum_terms[i-1]) < 1e-5)
      break
  }
  if(abs(sum_terms[N - 1] - sum_terms[N]) > 1e-5)
    warning("Series not yet converged at N = 250.")
  p + sum(sum_terms[1:i])
}

#' Density of the Kimura distribution
#'
#' Estimates the probability density function at x for heteroplasmy under no
#' selection pressure.
#'
#' @param x number in (0; 1): the point at which F() is calculated.
#' @inheritParams .f0
#'
#' @return The density of Kimura(p, d), evaluated at x.
#' @export
#'
#' @examples
#' .phi(0.2, 0.3, 0.5)
.phi <- function(x, p, b) {
  q <- 1 - p
  N <- 250
  sum_terms <- numeric(N)
  i_coef <- sapply(1:N, function(i) i * (i + 1) * (2 * i + 1))
  i_coef <- p * q * i_coef
  b_coef <- b^choose(2:(N + 1), 2)
  x_hg <- .hypgeo(1, x)
  p_hg <- .hypgeo(1, p)
  sum_terms[1] <- p_hg * x_hg * i_coef[1] * b_coef[1]
  for (i in 2:N) {
    p_hg <- .hypgeo(i, p)
    x_hg <- .hypgeo(i, x)
    sum_terms[i] <- p_hg * x_hg * i_coef[i] * b_coef[i]
    if (abs(sum_terms[i]-sum_terms[i-1]) < 1e-5)
      break
  }
  if(abs(sum_terms[N - 1] - sum_terms[N]) > 1e-5)
    warning("Series not yet converged at N = 250.")
  sum(sum_terms[1:i])
}

#' The Kimura distribution
#'
#' Estimates the Kimura probability density function.
#'
#' @inheritParams .phi
#'
#' @return The density of Kimura(p, d), evaluated at x.
#' @export
#'
#' @examples
#' dkimura(0, 0.3, 0.5)
#' dkimura(0.2, 0.3, 0.5)
#' dkimura(1, 0.3, 0.5)
dkimura <- function(x, p, b) {
  if (x == 0)
    return(.f0(p, b))
  if (x == 1)
    return(.f1(p, b))
  .phi(x, p, b)
}

#' The Kimura distribution
#'
#' Estimates the Kimura cumulative distribution function.
#'
#' @inheritParams .phi
#'
#' @return The CDF of Kimura(p, d), evaluated at x.
#' @export
#'
#' @examples
#' pkimura(0, 0.3, 0.5)
#' pkimura(0.2, 0.3, 0.5)
#' pkimura(1, 0.3, 0.5)
pkimura <- function(x, p, b) {
  if (x == 0)
    return(.f0(p, b))
  if (x == 1)
    return(1)
  prob_loss <- .f0(p, b)
  dst <- seq(1e-3, x, 1e-3)
  if(dst[length(dst)] < x)
    dst <- c(dst, x)
  d <- length(dst)
  delta <- dst[d] - dst[d - 1]
  dst_y <- sapply(dst, function(x) .phi(x, p, b))
  prob_loss +
    1e-3 * .5 * (sum(dst_y[1:(d-1)]) + sum(dst_y[2:(d-2)])) +
    delta * .5 * (dst_y[d-1] + dst_y[d])
}

#' The Kimura distribution
#'
#' Generates from the Kimura distribution.
#'
#' @param n number of observations to generate.
#' @inheritParams .phi
#'
#' @return A vector of length n.
#' @export
#'
#' @examples
#' rkimura(10, 0.3, 0.5)
rkimura <- function(n, p, b) {
  dst_x <- seq(0, 1, 1e-3)
  dst_y <- sapply(dst_x, function(x) dkimura(x, p, b))
  spl <- sample(c(0, 0.5, 1),
                n,
                replace = TRUE,
                prob = c(.f0(p, b), 1 - .f0(p, b) - .f1(p, b), .f1(p, b)))
  idx <- which(spl == 0.5)
  m <- length(idx)
  spl[idx] <- stats::approx(dst_y, dst_x, stats::runif(m))$y
  spl
}
