#' Compute the moments of the relative range
#'
#' The relative range is `W = R / sigma`, the range of `n` independent standard normal values.
#' Its upper tail is `P(W > w) = 1 - ptukey(w, nmeans = n, df = Inf)`, so
#' `E(W)` and `E(W^2)` are integrals of that tail.
#'
#' @param n Subgroup size.
#' @return Named numeric vector `c(m1 = E(W), m2 = E(W^2))`.
#' @noRd
compute_rel_range_moments <- function(n) {
  tail <- function(w) 1 - stats::ptukey(w, nmeans = n, df = Inf)
  m1 <- stats::integrate(tail, 0, Inf, rel.tol = 1e-10)$value
  m2 <- stats::integrate(function(w) 2 * w * tail(w), 0, Inf, rel.tol = 1e-10)$value
  c(m1 = m1, m2 = m2)
}

# Moments already computed in this R session, stored by subgroup size. The values depend only on
# `n`, so `d2()`, `d3()`, and `qc_constants()` share one computation per size.
.moment_cache <- new.env(parent = emptyenv())

#' Moments of the relative range, computed once per subgroup size
#'
#' @param n A single subgroup size (whole number of at least 2).
#' @return Named numeric vector `c(m1 = E(W), m2 = E(W^2))`.
#' @noRd
rel_range_moments <- function(n) {
  key <- as.character(n)
  moments <- .moment_cache[[key]]
  if (is.null(moments)) {
    moments <- compute_rel_range_moments(n)
    assign(key, moments, envir = .moment_cache)
  }
  moments
}

#' Check a subgroup size
#'
#' @param n Subgroup size(s): whole numbers of at least 2.
#' @return `n`, invisibly; an error otherwise.
#' @noRd
check_n <- function(n) {
  if (!is.numeric(n) || anyNA(n) || any(n < 2) || any(n != round(n))) {
    stop("`n` must be whole numbers of at least 2 (the subgroup size).", call. = FALSE)
  }
  invisible(n)
}

#' The constant d2: mean of the relative range
#'
#' Let \eqn{Z_1, \ldots, Z_n}{Z_1, ..., Z_n} be independent standard normal values and let
#' \deqn{W = \max_j Z_j - \min_j Z_j}{W = max_j Z_j - min_j Z_j}
#' be their range. For a sample of size \eqn{n} from a normal process with standard deviation
#' \eqn{\sigma}{sigma}, the sample range is \eqn{R = \sigma W}{R = sigma W}, so \eqn{W} is called the *relative
#' range*. The constant \eqn{d_2} is its mean:
#' \deqn{d_2 = E(W) = \int_0^\infty P(W > w)\, dw, \qquad P(W \le w) = \texttt{ptukey}(w, n, \infty).}{d_2 = E(W) = integral_0^Inf P(W > w) dw, P(W <= w) = ptukey(w, n, Inf).}
#' Therefore
#' \deqn{E(R) = d_2\,\sigma,}{E(R) = d_2 sigma,}
#' which gives the estimator of \eqn{\sigma}{sigma} from the average subgroup range \eqn{\bar R}{Rbar}:
#' \deqn{\hat\sigma = \frac{\bar R}{d_2}.}{sigma-hat = (Rbar)/(d_2).}
#' The constant also gives the \eqn{\bar X}{Xbar} chart factor \eqn{A_2 = 3/(d_2\sqrt{n})}{A_2 = 3/(d_2*sqrt(n))}, so that the
#' limits are \eqn{\bar{\bar x} \pm A_2 \bar R}{xbarbar +/- A_2 Rbar}.
#'
#' The value is computed by numerical integration, so it is exact for any subgroup size
#' \eqn{n \ge 2}{n >= 2} (a published table is rounded, usually to 3 decimals). For \eqn{n = 2},
#' \eqn{d_2 = 2/\sqrt{\pi}}{d_2 = 2/sqrt(pi)}. It increases with \eqn{n}: a larger sample has a wider expected range.
#'
#' @param n Subgroup size: whole numbers of at least 2 (a vector is allowed).
#' @return A numeric vector the same length as `n`.
#' @family control chart constants
#' @seealso [qc_constants()] for all constants in one table, and [qc_chart()], which uses
#'   \eqn{\bar R/d_2}{Rbar/d_2} to estimate \eqn{\sigma}{sigma}.
#' @examples
#' d2(4)
#' d2(2:6)
#' d2(2) - 2 / sqrt(pi)    # exact value for n = 2
#'
#' # d2 is the mean range of n standard normal values: check by simulation
#' set.seed(1)
#' mean(replicate(10000, diff(range(rnorm(4)))))
#'
#' # Estimating sigma from an average range of 3.82 with subgroups of size 4
#' 3.82 / d2(4)
#' @export
d2 <- function(n) {
  check_n(n)
  vapply(n, function(k) unname(rel_range_moments(k)["m1"]), numeric(1))
}

#' The constant d3: standard deviation of the relative range
#'
#' With \eqn{W} the range of \eqn{n} independent standard normal values (see [d2()]), the
#' constant \eqn{d_3} is the standard deviation of \eqn{W}:
#' \deqn{d_3 = \sqrt{E(W^2) - d_2^2}, \qquad E(W^2) = \int_0^\infty 2w\, P(W > w)\, dw.}{d_3 = sqrt(E(W^2) - d_2^2), E(W^2) = integral_0^Inf 2w P(W > w) dw.}
#' Since \eqn{R = \sigma W}{R = sigma W},
#' \deqn{SD(R) = d_3\,\sigma.}{SD(R) = d_3 sigma.}
#' A three-sigma control chart for the range therefore has center line \eqn{\bar R = d_2\hat\sigma}{Rbar = d_2*sigma-hat}
#' and limits
#' \deqn{\bar R \pm 3 d_3\hat\sigma = \bar R\left(1 \pm \frac{3d_3}{d_2}\right),}{Rbar +/- 3 d_3*sigma-hat = Rbar(1 +/- (3d_3)/(d_2)),}
#' which defines the constants
#' \deqn{D_4 = 1 + \frac{3d_3}{d_2}, \qquad D_3 = \max\left(0,\, 1 - \frac{3d_3}{d_2}\right).}{D_4 = 1 + (3d_3)/(d_2), D_3 = max(0, 1 - (3d_3)/(d_2)).}
#' For \eqn{n \le 6}{n <= 6} the lower limit is negative and is set to 0, so \eqn{D_3 = 0}.
#'
#' For \eqn{n = 2}, \eqn{d_3 = \sqrt{2(1 - 2/\pi)}}{d_3 = sqrt(2(1 - 2/pi))}. The values are computed by numerical
#' integration, so any subgroup size \eqn{n \ge 2}{n >= 2} is allowed.
#'
#' @inheritParams d2
#' @return A numeric vector the same length as `n`.
#' @family control chart constants
#' @seealso [qc_constants()] for all constants in one table, including \eqn{D_3} and \eqn{D_4}.
#' @examples
#' d3(4)
#' d3(2:6)
#' d3(2) - sqrt(2 * (1 - 2 / pi))    # exact value for n = 2
#'
#' # d3 is the SD of the range of n standard normal values: check by simulation
#' set.seed(1)
#' sd(replicate(10000, diff(range(rnorm(4)))))
#'
#' # R chart factors for n = 5, built from d2 and d3
#' c(D3 = max(0, 1 - 3 * d3(5) / d2(5)), D4 = 1 + 3 * d3(5) / d2(5))
#' @export
d3 <- function(n) {
  check_n(n)
  vapply(n, function(k) {
    m <- rel_range_moments(k)
    sqrt(m[["m2"]] - m[["m1"]]^2)
  }, numeric(1))
}

#' The constant c4: bias correction for the sample standard deviation
#'
#' For a sample of size \eqn{n} from a normal process with standard deviation \eqn{\sigma}{sigma}, the
#' sample standard deviation \eqn{s} is biased: \eqn{E(s) = c_4\,\sigma}{E(s) = c_4 sigma}, where
#' \deqn{c_4 = \sqrt{\frac{2}{n - 1}}\;\frac{\Gamma(n/2)}{\Gamma((n-1)/2)}.}{c_4 = sqrt((2)/(n - 1)) (Gamma(n/2))/(Gamma((n-1)/2)).}
#' (The sample *variance* \eqn{s^2} is unbiased for \eqn{\sigma^2}{sigma^2}, but taking the square root
#' makes \eqn{s} biased low, so \eqn{c_4 < 1}.) Given the average subgroup SD \eqn{\bar s}{sbar}, an
#' unbiased estimator of \eqn{\sigma}{sigma} is
#' \deqn{\hat\sigma = \frac{\bar s}{c_4}.}{sigma-hat = (sbar)/(c_4).}
#' The constant also gives the \eqn{\bar X}{Xbar} chart factor \eqn{A_3 = 3/(c_4\sqrt{n})}{A_3 = 3/(c_4*sqrt(n))} for charts
#' based on \eqn{\bar s}{sbar}, so that the limits are \eqn{\bar{\bar x} \pm A_3\bar s}{xbarbar +/- A_3*sbar}.
#'
#' For \eqn{n = 2}, \eqn{c_4 = \sqrt{2/\pi}}{c_4 = sqrt(2/pi)}. The constant increases toward 1 as \eqn{n} grows.
#'
#' @inheritParams d2
#' @return A numeric vector the same length as `n`.
#' @family control chart constants
#' @seealso [qc_constants()] for all constants in one table, including \eqn{A_3}, \eqn{B_3},
#'   and \eqn{B_4}.
#' @examples
#' c4(5)
#' c4(c(2, 5, 10, 25))
#' c4(2) - sqrt(2 / pi)    # exact value for n = 2
#'
#' # c4 is the expected sample SD when sigma = 1: check by simulation
#' set.seed(1)
#' mean(replicate(10000, sd(rnorm(5))))
#'
#' # Unbiased estimate of sigma from an average subgroup SD of 1.6 with n = 5
#' 1.6 / c4(5)
#' @export
c4 <- function(n) {
  check_n(n)
  sqrt(2 / (n - 1)) * exp(lgamma(n / 2) - lgamma((n - 1) / 2))
}

#' Table of control chart constants
#'
#' One row per subgroup size, with all the constants used for the \eqn{\bar X}{Xbar}, \eqn{R}, and
#' \eqn{S} chart limits. They are built from [d2()], [d3()], and [c4()]:
#' \deqn{A_2 = \frac{3}{d_2\sqrt n}, \quad D_3 = \max\left(0,\, 1 - \frac{3d_3}{d_2}\right), \quad
#' D_4 = 1 + \frac{3d_3}{d_2},}{A_2 = (3)/(d_2*sqrt(n)), D_3 = max(0, 1 - (3d_3)/(d_2)), D_4 = 1 + (3d_3)/(d_2),}
#' \deqn{A_3 = \frac{3}{c_4\sqrt n}, \quad B_3 = \max\left(0,\, 1 - \frac{3c_5}{c_4}\right), \quad
#' B_4 = 1 + \frac{3c_5}{c_4},}{A_3 = (3)/(c_4*sqrt(n)), B_3 = max(0, 1 - (3c_5)/(c_4)), B_4 = 1 + (3c_5)/(c_4),}
#' where \eqn{c_5 = \sqrt{1 - c_4^2}}{c_5 = sqrt(1 - c_4^2)} is the standard deviation of \eqn{s/\sigma}{s/sigma}.
#'
#' Which constants go with which chart (center line CL, limits LCL and UCL):
#'
#' | Chart | Estimate of \eqn{\sigma}{sigma} | CL | LCL, UCL |
#' |:--|:--|:--|:--|
#' | \eqn{\bar X}{Xbar} (with \eqn{R}) | \eqn{\bar R/d_2}{Rbar/d_2} | \eqn{\bar{\bar x}}{xbarbar} | \eqn{\bar{\bar x} \mp A_2\bar R}{xbarbar -/+ A_2*Rbar} |
#' | \eqn{R} | \eqn{\bar R/d_2}{Rbar/d_2} | \eqn{\bar R}{Rbar} | \eqn{D_3\bar R,\; D_4\bar R}{D_3*Rbar, D_4*Rbar} |
#' | \eqn{\bar X}{Xbar} (with \eqn{S}) | \eqn{\bar s/c_4}{sbar/c_4} | \eqn{\bar{\bar x}}{xbarbar} | \eqn{\bar{\bar x} \mp A_3\bar s}{xbarbar -/+ A_3*sbar} |
#' | \eqn{S} | \eqn{\bar s/c_4}{sbar/c_4} | \eqn{\bar s}{sbar} | \eqn{B_3\bar s,\; B_4\bar s}{B_3*sbar, B_4*sbar} |
#'
#' @param n Subgroup sizes (whole numbers of at least 2).
#' @return A tibble with columns `n`, `d2`, `d3`, `c4`, `A2`, `A3`, `B3`, `B4`, `D3`, `D4`.
#' @family control chart constants
#' @examples
#' qc_constants(2:6)
#'
#' # limits of an X-bar chart with xbarbar = 50 and Rbar = 2 for subgroups of size 5
#' 50 + c(-1, 1) * qc_constants(5)$A2 * 2
#' @export
qc_constants <- function(n = 2:10) {
  check_n(n)
  d2n <- d2(n)
  d3n <- d3(n)
  c4n <- c4(n)
  c5n <- sqrt(1 - c4n^2)
  tibble::tibble(
    n = n, d2 = d2n, d3 = d3n, c4 = c4n,
    A2 = 3 / (d2n * sqrt(n)),
    A3 = 3 / (c4n * sqrt(n)),
    B3 = pmax(0, 1 - 3 * c5n / c4n),
    B4 = 1 + 3 * c5n / c4n,
    D3 = pmax(0, 1 - 3 * d3n / d2n),
    D4 = 1 + 3 * d3n / d2n
  )
}
