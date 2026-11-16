# Gauss-Legendre nodes and weights on [-1, 1] (Golub-Welsch).
gauss_legendre <- function(n) {
  i <- seq_len(n - 1)
  b <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1)] <- b
  J[cbind(i + 1, i)] <- b
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = 2 * e$vectors[1, ]^2)
}

check_shift <- function(shift) {
  if (!is.numeric(shift) || anyNA(shift) || any(!is.finite(shift))) stop("`shift` must be numeric with no missing values.", call. = FALSE)
}

#' Average run length of Shewhart, CUSUM, and EWMA charts
#'
#' Exact (numerical) average run lengths of the charts of this package for **normal, independent** results with a shift
#' of the mean of `shift` standard deviations (`shift = 0` is the in-control ARL). All of them assume that the target
#' and the standard deviation are known.
#' * `shewhart_arl()`: the ARL of a Shewhart chart with limits at \eqn{\pm L} standard deviations of the plotted values,
#'   \eqn{1/p} with \eqn{p = P(|Z| > L)}, a geometric run length.
#' * `cusum_arl()`: the ARL of the tabular CUSUM of [cusum_chart()] with reference value `k` and decision interval `h`.
#'   Each one-sided CUSUM is solved from the integral equation of Goel and Wu (1971) by Gauss-Legendre quadrature with `n` nodes;
#'   for two sides the usual approximation \eqn{1/\text{ARL} = 1/\text{ARL}^+ + 1/\text{ARL}^-}{1/ARL = 1/ARL+ + 1/ARL-} is used,
#'   which is accurate to a few tenths of a percent.
#' * `ewma_arl()`: the ARL of the EWMA chart of [ewma_chart()] with the **asymptotic** limits, from the integral equation of
#'   Crowder (1987) by Gauss-Legendre quadrature with `n` nodes.
#'
#' @param shift Shift of the mean in standard deviations (a number or a vector).
#' @param L Width of the limits (Shewhart: in standard deviations of the plotted values; EWMA: in standard deviations of the EWMA).
#' @param k,h Reference value and decision interval of the CUSUM, in standard deviations.
#' @param sides `2` (default) for the two-sided CUSUM, or `1` for the upper one-sided CUSUM.
#' @param head_start Starting value of the CUSUM (default 0).
#' @param lambda Weight of the EWMA, in (0, 1].
#' @param n Number of quadrature nodes. Default 60.
#' @return A numeric vector of average run lengths, one per `shift`.
#' @references Goel, A. L. and Wu, S. M. (1971). Determination of A.R.L. and a contour nomogram for CUSUM charts to control normal mean.
#'   *Technometrics*, 13, 221-230. Crowder, S. V. (1987). A simple method for studying run-length distributions of exponentially
#'   weighted moving average charts. *Technometrics*, 29, 401-407.
#' @seealso [find_cusum_h()], [find_ewma_L()], [sim_arl()]
#' @examples
#' shewhart_arl(c(0, 0.5, 1, 2))
#' cusum_arl(c(0, 0.5, 1, 2), k = 0.5, h = 5)      # in-control ARL 465
#' ewma_arl(c(0, 0.5, 1, 2), lambda = 0.2, L = 2.86)
#' @export
shewhart_arl <- function(shift = 0, L = 3) {
  check_shift(shift); check_number(L, "L", positive = TRUE)
  p <- stats::pnorm(-L - shift) + stats::pnorm(-L + shift)
  1 / p
}

cusum_arl_one <- function(shift, k, h, n, head_start) {
  gl <- gauss_legendre(n)
  x <- c(0, h * (1 + gl$x) / 2)
  w <- c(0, h * gl$w / 2)                       # the point 0 only carries the head-start-free equation
  A <- stats::pnorm(k - x - shift)               # P(next sum is 0)
  K <- outer(x, x[-1], function(u, y) stats::dnorm(y - u + k - shift))
  K <- sweep(K, 2, w[-1], `*`)
  M <- diag(n + 1) - cbind(A, K)
  v <- tryCatch(solve(M, rep(1, n + 1)), error = function(e) NULL)
  if (is.null(v) || !is.finite(v[1]) || v[1] > 1e9 || v[1] < 1) return(Inf)    # a shift away from the limit: the ARL is astronomically large
  if (head_start == 0) return(v[1])
  u <- head_start
  1 + v[1] * stats::pnorm(k - u - shift) + sum(w[-1] * stats::dnorm(x[-1] - u + k - shift) * v[-1])
}

#' @rdname shewhart_arl
#' @export
cusum_arl <- function(shift = 0, k = 0.5, h = 5, sides = 2, head_start = 0, n = 60) {
  check_shift(shift); check_number(k, "k", positive = TRUE); check_number(h, "h", positive = TRUE)
  if (!sides %in% c(1, 2)) stop("`sides` must be 1 or 2.", call. = FALSE)
  check_finite(head_start, "head_start")
  if (head_start < 0 || head_start >= h) stop("`head_start` must be between 0 and `h`.", call. = FALSE)
  check_count(n, "n", minimum = 10)
  vapply(shift, function(d) {
    up <- cusum_arl_one(d, k, h, n, head_start)
    if (sides == 1) up else 1 / (1 / up + 1 / cusum_arl_one(-d, k, h, n, head_start))
  }, numeric(1))
}

#' @rdname shewhart_arl
#' @export
ewma_arl <- function(shift = 0, lambda = 0.2, L = 3, n = 60) {
  check_shift(shift); check_number(L, "L", positive = TRUE); check_count(n, "n", minimum = 10)
  if (!(is.numeric(lambda) && length(lambda) == 1 && is.finite(lambda) && lambda > 0 && lambda <= 1)) {
    stop("`lambda` must be a single number in (0, 1].", call. = FALSE)
  }
  cc <- L * sqrt(lambda / (2 - lambda))
  gl <- gauss_legendre(n)
  x <- cc * gl$x; w <- cc * gl$w
  vapply(shift, function(d) {
    K <- outer(x, x, function(u, y) stats::dnorm((y - (1 - lambda) * u) / lambda - d) / lambda)
    K <- sweep(K, 2, w, `*`)
    v <- tryCatch(solve(diag(n) - K, rep(1, n)), error = function(e) NULL)
    if (is.null(v)) return(Inf)
    arl <- 1 + sum(w * stats::dnorm(x / lambda - d) / lambda * v)
    if (!is.finite(arl) || arl > 1e9 || arl < 1) Inf else arl
  }, numeric(1))
}

#' Find the decision interval or the limit width for a target in-control ARL
#'
#' `find_cusum_h()` finds the decision interval `h` of the CUSUM with reference value `k`, and `find_ewma_L()` the width `L`
#' of the EWMA limits with weight `lambda`, that give the in-control average run length `arl0` (see [cusum_arl()] and
#' [ewma_arl()]). Matching the in-control ARL of different charts is how their detection speeds are compared fairly.
#'
#' @param arl0 Wanted in-control average run length (greater than 2).
#' @param k Reference value of the CUSUM. Default 0.5.
#' @param sides `2` (default) or `1`.
#' @param lambda Weight of the EWMA. Default 0.2.
#' @return A single number.
#' @seealso [cusum_arl()], [ewma_arl()]
#' @examples
#' find_cusum_h(370, k = 0.5)
#' find_ewma_L(370, lambda = 0.2)
#' @export
find_cusum_h <- function(arl0, k = 0.5, sides = 2) {
  check_number(arl0, "arl0", positive = TRUE); check_number(k, "k", positive = TRUE)
  if (arl0 <= 2) stop("`arl0` must be greater than 2.", call. = FALSE)
  f <- function(h) min(cusum_arl(0, k, h, sides), 1e9) - arl0
  upper <- 25
  if (f(upper) < 0) stop("The decision interval for this ARL is above 25.", call. = FALSE)
  stats::uniroot(f, c(0.05, upper), tol = 1e-8)$root
}

#' @rdname find_cusum_h
#' @export
find_ewma_L <- function(arl0, lambda = 0.2) {
  check_number(arl0, "arl0", positive = TRUE)
  if (arl0 <= 2) stop("`arl0` must be greater than 2.", call. = FALSE)
  f <- function(L) min(ewma_arl(0, lambda, L), 1e9) - arl0
  if (f(6) < 0) stop("The limit width for this ARL is above 6.", call. = FALSE)
  stats::uniroot(f, c(0.2, 6), tol = 1e-8)$root
}
