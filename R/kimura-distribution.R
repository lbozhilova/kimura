#' Calculate \eqn{F(i-1, i+2, 2, x)}
#'
#' Use Gauss' contiguous relations to calculate the hypergeometric function
#' \eqn{F(i-1, i+2, 2, x)} needed for estimating the Kimura distribution.
#'
#' @param i positive integer.
#' @param x number in \eqn{(0; 1)}: the point at which the function is
#'   calculated.
#'
#' @return The value of \eqn{F(i-1, i+2, 2, x)}, calculated recursively.
#'
#' @references Wonnapinij, Passorn, Patrick F. Chinnery, and David C. Samuels.
#'   "The distribution of mitochondrial DNA heteroplasmy due to random genetic
#'   drift." The American Journal of Human Genetics 83.5 (2008): 582-593.
#'
#'   Hoang-Binh, D. "A program to compute exact hydrogenic radial integrals,
#'   oscillator strengths, and Einstein coefficients, for principal quantum
#'   numbers up to n≈ 1000." Computer Physics Communications 166.3 (2005):
#'   191-196.
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
#' Estimate the probability that an allele is lost due to random genetic drift.
#'
#' @param p number in \eqn{(0; 1)}: initial heteroplasmy.
#' @param b number in \eqn{(0; 1)}: drift parameter.
#'
#' @details The parameter \code{p} corresponds to the initial heteroplasmy, and
#'   therefore is also equal to the mean heteroplasmy after genetic drift.
#'
#'   The parameter \code{b} corresponds to the amount of drift, and is equal to
#'   \code{b = exp(-t/N_eff)}, where \code{t} denotes time and \code{N_eff}
#'   denotes effective sample size.
#'
#' @return The probability the allele is lost due to genetic drift.
#'
#' @references Wonnapinij, Passorn, Patrick F. Chinnery, and David C. Samuels.
#'   "The distribution of mitochondrial DNA heteroplasmy due to random genetic
#'   drift." The American Journal of Human Genetics 83.5 (2008): 582-593.
#'
#'   Kimura, Motoo. "Solution of a process of random genetic drift with a
#'   continuous model." Proceedings of the National Academy of Sciences of the
#'   United States of America 41.3 (1955): 144.
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
    if (abs(sum_terms[i]-sum_terms[i-1]) < 1e-4)
      break
  }
  if(abs(sum_terms[N - 1] - sum_terms[N]) > 1e-4)
    warning("Series not yet converged at N = 250.")
  max(q + sum(sum_terms[1:i]), 0)
}

#' Probability of fixing allele
#'
#' Estimate the probability that an allele is fixed (i.e. the wildtype is lost)
#' due to random genetic drift.
#'
#' @inheritParams .f0
#'
#' @details The parameter \code{p} corresponds to the initial heteroplasmy, and
#'   therefore is also equal to the mean heteroplasmy after genetic drift.
#'
#'   The parameter \code{b} corresponds to the amount of drift, and is equal to
#'   \code{b = exp(-t/N_eff)}, where \code{t} denotes time and \code{N_eff}
#'   denotes effective sample size.
#'
#' @return The probability the allele is fixed due to genetic drift.
#'
#' @references Wonnapinij, Passorn, Patrick F. Chinnery, and David C. Samuels.
#'   "The distribution of mitochondrial DNA heteroplasmy due to random genetic
#'   drift." The American Journal of Human Genetics 83.5 (2008): 582-593.
#'
#'   Kimura, Motoo. "Solution of a process of random genetic drift with a
#'   continuous model." Proceedings of the National Academy of Sciences of the
#'   United States of America 41.3 (1955): 144.
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
    if (abs(sum_terms[i]-sum_terms[i-1]) < 1e-4)
      break
  }
  if(abs(sum_terms[N - 1] - sum_terms[N]) > 1e-4)
    warning("Series not yet converged at N = 250.")
  max(p + sum(sum_terms[1:i]), 0)
}

#' Density of the Kimura distribution
#'
#' Estimate the probability density function at \code{x} for heteroplasmy under
#' no selection pressure.
#'
#' @param x number in \eqn{(0; 1)}: the point at which the density is
#'   calculated.
#' @inheritParams .f0
#'
#' @return The density of \eqn{Kimura(p, b)}, evaluated at \eqn{x}.
#'
#' @references Wonnapinij, Passorn, Patrick F. Chinnery, and David C. Samuels.
#'   "The distribution of mitochondrial DNA heteroplasmy due to random genetic
#'   drift." The American Journal of Human Genetics 83.5 (2008): 582-593.
#'
#'   Kimura, Motoo. "Solution of a process of random genetic drift with a
#'   continuous model." Proceedings of the National Academy of Sciences of the
#'   United States of America 41.3 (1955): 144.
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
    if (abs(sum_terms[i]-sum_terms[i-1]) < 1e-4)
      break
  }
  if(abs(sum_terms[N - 1] - sum_terms[N]) > 1e-4)
    warning("Series not yet converged at N = 250.")
  max(sum(sum_terms[1:i]), 0)
}

#' The Kimura distribution
#'
#' Estimate the \eqn{Kimura(p, b)} probability density function.
#'
#' @param x numeric vector of values in \eqn{[0; 1]}.
#' @inheritParams .phi
#'
#' @return The density of \eqn{Kimura(p, b)}, evaluated at \code{x}. Note there
#'   is a fixed probability at both extremes of the distribution.
#'
#' @references Wonnapinij, Passorn, Patrick F. Chinnery, and David C. Samuels.
#'   "The distribution of mitochondrial DNA heteroplasmy due to random genetic
#'   drift." The American Journal of Human Genetics 83.5 (2008): 582-593.
#'
#'   Kimura, Motoo. "Solution of a process of random genetic drift with a
#'   continuous model." Proceedings of the National Academy of Sciences of the
#'   United States of America 41.3 (1955): 144.
#' @export
#'
#' @examples
#' dkimura(c(0, 0.2, 0.3), 0.3, 0.5)
dkimura <- function(x, p, b) {
  if (b == 0)
    return(stats::dbinom(x, 1, p))
  if (b == 1)
    return(1 * x == p)
  y <- numeric(length(x))
  y[x == 0] <- .f0(p, b)
  y[x == 1] <- .f1(p, b)
  z <- x[x > 0 & x < 1]
  lz <- length(z)
  if (lz == 0)
    return(y)
  q <- 1 - p
  N <- 250
  sum_terms <- matrix(0, ncol = lz, nrow = N)
  i_coef <- sapply(1:N, function(i) i * (i + 1) * (2 * i + 1))
  i_coef <- p * q * i_coef
  b_coef <- b^choose(2:(N + 1), 2)
  p_hg <- .hypgeo(1, p)
  z_hg <- sapply(z, function(x) .hypgeo(1, x))
  sum_terms[1, ] <- p_hg * z_hg * i_coef[1] * b_coef[1]
  j <- 1:lz
  for (i in 2:N) {
    z_hg <- sapply(z[j], function(x) .hypgeo(i, x))
    p_hg <- .hypgeo(i, p)
    sum_terms[i, j] <- p_hg * z_hg * i_coef[i] * b_coef[i]
    done_idx <- abs(sum_terms[i-1, j] - sum_terms[i, j]) < 1e-4
    j <- j[!done_idx]
    if (length(j) == 0)
      break
  }
  if(any(abs(sum_terms[N - 1, ] - sum_terms[N, ]) > 1e-4))
    warning("Series not yet converged at N = 250.")
  y[x > 0 & x < 1] <- colSums(sum_terms)
  y[y < 0] <- 0
  y
}

#' The Kimura distribution
#'
#' Estimate the \eqn{Kimura(p, b)} cumulative distribution function.
#'
#' @inheritParams dkimura
#'
#' @return The cumulative distribution function of \eqn{Kimura(p, b)}, evaluated
#'   at \code{x}. Note there is a fixed probability at both extremes of the
#'   distribution. the distribution is estimated by evaluating \code{dkimura(p,
#'   b)} at regular intervals of length \code{1e-4} and applying the trapezoidal
#'   rule.
#'
#' @references Wonnapinij, Passorn, Patrick F. Chinnery, and David C. Samuels.
#'   "The distribution of mitochondrial DNA heteroplasmy due to random genetic
#'   drift." The American Journal of Human Genetics 83.5 (2008): 582-593.
#' @export
#'
#' @examples
#' pkimura(c(0, 0.2, 0.3), 0.3, 0.5)
pkimura <- function(x, p, b) {
  if (b == 0)
    return(stats::pbinom(x, 1, p))
  if (b == 1)
    return(1 * (x >= p))
  y <- numeric(length(x))
  y[x < 0] <- 0
  y[x == 0] <- .f0(p, b)
  y[x >= 1] <- 1
  z <- x[x > 0 & x < 1]
  lz <- length(z)
  if (lz == 0)
    return(y)
  density_z <- dkimura(z, p, b)
  k <- seq(1e-4, max(z) + 1e-4, 1e-4)
  density_k <- dkimura(k, p, b)
  cdf_k <- .pkimura_full(p, b)[-c(1, 10001)]
  get_cdf <- function(ix) {
    x <- z[ix]
    ik <- floor(x * 1e4)
    if(ik == 0)
      return(.f0(p, b) + 0.5 * x * density_z[ix])
    if((ik+1) == 1e4) {
      d <- 1 - x
      return(1 - .f1(p, b) - 0.5 * d * density_z[ix])
    }
    d <- x - k[ik]
    cdf_k[ik] + 0.5 * d * (density_k[ik] + density_z[ix])
  }
  y[x > 0 & x < 1] <- sapply(1:length(z), get_cdf)
  y
}

#' The Kimura distribution
#'
#' Generate random samples from the \eqn{Kimura(p, b)} distribution.
#'
#' @param n number of observations to generate.
#' @inheritParams .phi
#'
#' @return A vector of length \code{n} of heteroplasmy fractions generated from
#'   the \eqn{Kimura(p, b)} distribution. Note that while most real-life
#'   heteroplasmy fraction data is rounded to 2 significant digits, these
#'   randomly generated values are not.
#'
#' @references Wonnapinij, Passorn, Patrick F. Chinnery, and David C. Samuels.
#'   "The distribution of mitochondrial DNA heteroplasmy due to random genetic
#'   drift." The American Journal of Human Genetics 83.5 (2008): 582-593.
#' @export
#'
#' @examples
#' rkimura(10, 0.3, 0.5)
rkimura <- function(n, p, b) {
  if (b == 0)
    return(stats::rbinom(n, 1, p))
  if (b == 1)
    return(rep(p, n))
  x <- seq(0, 1, 1e-4)
  cdf_x <- .pkimura_full(p, b)
  spl <- sample(c(0, 0.5, 1), n, replace = TRUE,
                prob = c(.f0(p, b), 1 - .f0(p, b) - .f1(p, b), .f1(p, b)))
  idx <- which(spl == 0.5)
  m <- length(idx)
  spl[idx] <- stats::approx(cdf_x, x, stats::runif(m, min = cdf_x[1]))$y
  spl
}
