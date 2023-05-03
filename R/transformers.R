#' Normalised heteroplasmy variance
#'
#' @param H numeric vector of heteroplasmy fractions, either in \eqn{[0; 1]} or
#'    represented as percentages in \eqn{[0; 100]}
#' @param percent logical; are values in H represented as percentage points?
#'
#' @return The normalised heteroplasmy variance, which is typically in
#'    \eqn{[0; 1]}. Lower values imply lower variance and less genetic
#'    drift, with 0 corresponding to all values of H being equal.
#'
#'    Note \code{nvar(H)} is equivalent to \code{1-estimate_parameters(H)["b"]}.
#' @export
#'
#' @examples
#' # Using percentage points
#' H <- sample(0:100, 10)
#' nvar(H)
#'
#' # Using heteroplasmy fractions between 0 and 1
#' h <- stats::runif(10)
#' nvar(h, percent = FALSE)
nvar <- function(H, percent = TRUE) {
  if (percent) H <- .01 * H
  stats::var(H) / (mean(H) * (1 - mean(H)))
}

#' Normalised heteroplasmy shift
#'
#' @inheritParams nvar
#' @param H0 numeric vector of baseline heteroplasmy fractions, either in
#'    \eqn{[0; 1]} or represented as percentages in \eqn{[0; 100]}
#' @param tol tolerance with which to resolve homoplasmy
#'
#' @return A vector of normalised heteroplasmy shift calculated as
#'
#' \eqn{shift(H; H_0) = \frac{H(1-H_0)}{H_0(1-H)}}.
#'
#' Homoplasmic states are resolved as \code{tol} for homoplasmic wildtype, and
#' \code{(1 - tol)} for homoplasmic mutant.
#'
#' @export
#'
#' @examples
#' # Using percentage points
#' H <- sample(0:100, 10)
#' H0 <- 50
#' shift(H, H0)
#'
#' # Using fractions
#' h <- stats::runif(10)
#' h0 <- 0.3
#' shift(h, h0, percent = FALSE)
shift <- function(H, H0, percent = TRUE, tol = 1e-4) {
  if (length(H0) != 1 & length(H0) != length(H))
    stop("The baseline heteropalsmy vector H0 should be either of length 1 or of the same length as H.")
  if (percent) {
    H <- .01 * H
    H0 <- .01 * H0
  }
  H[H == 0] <- tol
  H[H == 1] <- 1 - tol
  H0[H0 == 0] <- tol
  H0[H0 == 1] <- 1 - tol
  stats::qlogis(H) - stats::qlogis(H0)
}
