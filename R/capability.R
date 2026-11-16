#' Process capability indices with confidence intervals
#'
#' Computes the capability indices of a sample of individual measurements against specification limits,
#' with approximate confidence intervals, the observed and the normal-theory fraction outside the
#' specification, and the defects per million.
#'
#' With \eqn{\hat\sigma} the standard deviation used for the **capability** (by default the sample
#' standard deviation \eqn{s}; a within-subgroup estimate such as \eqn{\bar R/d_2}{Rbar/d2} can be supplied with `sigma`) and
#' \eqn{\bar x}{xbar} the sample mean,
#' * \eqn{C_p = (USL - LSL)/(6\hat\sigma)}{Cp = (USL - LSL)/(6 sigma-hat)}: the potential capability, which ignores centering;
#' * \eqn{C_{pl} = (\bar x - LSL)/(3\hat\sigma)}{Cpl = (xbar - LSL)/(3 sigma-hat)} and \eqn{C_{pu} = (USL - \bar x)/(3\hat\sigma)}{Cpu = (USL - xbar)/(3 sigma-hat)};
#' * \eqn{C_{pk} = \min(C_{pl}, C_{pu})}{Cpk = min(Cpl, Cpu)}: the capability including centering;
#' * \eqn{C_{pm} = (USL - LSL)/(6\sqrt{\hat\sigma^2 + (\bar x - T)^2})}{Cpm = (USL - LSL)/(6 sqrt(sigma-hat^2 + (xbar - T)^2))} when a `target` \eqn{T} is given.
#'
#' The **performance indices** \eqn{P_p} and \eqn{P_{pk}} are the same formulas with the overall standard
#' deviation \eqn{s} of all the results. They differ from \eqn{C_p} and \eqn{C_{pk}} only when `sigma` is supplied
#' (a within-subgroup estimate), in which case the difference measures the variation between subgroups.
#'
#' The confidence intervals assume normal, independent data. For \eqn{C_p} the interval is exact
#' (chi-square with `df` degrees of freedom). For \eqn{C_{pk}} (and \eqn{C_{pl}}, \eqn{C_{pu}}) `ci = "normal"` uses the
#' large-sample approximation of Bissell (1990),
#' \deqn{\hat C_{pk} \pm z\sqrt{\frac{1}{9n} + \frac{\hat C_{pk}^2}{2\,df}},}{Cpk-hat +/- z sqrt(1/(9n) + Cpk-hat^2/(2 df)),}
#' with \eqn{df = n - 1} when \eqn{\hat\sigma = s}; `ci = "bootstrap"` gives the percentile interval from `B` resamples
#' of the data (it needs `sigma = NULL`); `ci = "none"` skips the intervals. The one-sided lower
#' bound for \eqn{C_{pk}} (at the same confidence level) is also returned.
#'
#' @param x Numeric vector of individual measurements (no missing values, at least 3).
#' @param lsl,usl Lower and upper specification limits. Give one or both (\eqn{C_p} needs both).
#' @param target Optional target for \eqn{C_{pm}}.
#' @param sigma Optional standard deviation to use for the capability indices (for example a within-subgroup
#'   estimate). Default `NULL`: the sample standard deviation.
#' @param df Degrees of freedom of `sigma` for the intervals. Default `n - 1` when `sigma` is `NULL`; required
#'   (otherwise the intervals are `NA`) when `sigma` is supplied.
#' @param conf Confidence level. Default 0.95.
#' @param ci `"normal"` (default), `"bootstrap"`, or `"none"`.
#' @param B Number of bootstrap resamples. Default 2000. Use [set.seed()] for reproducibility.
#' @return An object of class `capability`: a list with `indices` (a tibble with the estimate and the lower and upper
#'   confidence limits of each index), `n`, `mean`, `sd`, `sigma` (the SD used), `lower_one_sided` (lower confidence
#'   bound for Cpk), `observed_out` (observed fraction outside the specification), `expected_out` (normal-theory fraction
#'   with the sample mean and the SD used for the indices, `sigma`), `dpmo` (expected defects per million, normal theory), `lsl`, `usl`, `conf`, and `ci`.
#' @references Bissell, A. F. (1990). How reliable is your capability index? *Applied Statistics*, 39, 331-340.
#' @seealso [sim_capability()], [dpmo()], [sigma_level()]
#' @examples
#' capability(tablet_supplier$weight_mg, lsl = 475, usl = 525)
#'
#' # Percentile bootstrap interval
#' set.seed(1)
#' capability(tablet_supplier$weight_mg, lsl = 475, usl = 525, ci = "bootstrap")
#'
#' # One-sided specification
#' capability(tablet_supplier$weight_mg, usl = 525)
#' @export
capability <- function(x, lsl = NULL, usl = NULL, target = NULL, sigma = NULL, df = NULL, conf = 0.95,
                       ci = c("normal", "bootstrap", "none"), B = 2000) {
  ci <- match.arg(ci)
  if (!is.numeric(x) || anyNA(x) || length(x) < 3) stop("`x` must be a numeric vector with at least 3 values and no missing values.", call. = FALSE)
  if (is.null(lsl) && is.null(usl)) stop("Give at least one specification limit (`lsl` or `usl`).", call. = FALSE)
  if (!is.null(lsl)) check_finite(lsl, "lsl")
  if (!is.null(usl)) check_finite(usl, "usl")
  if (!is.null(lsl) && !is.null(usl) && lsl >= usl) stop("`lsl` must be below `usl`.", call. = FALSE)
  if (!is.null(target)) check_finite(target, "target")
  check_probability(conf, "conf")
  n <- length(x); m <- mean(x); s <- stats::sd(x)
  if (s == 0) stop("The measurements have no variation.", call. = FALSE)
  if (!is.null(sigma)) {
    check_number(sigma, "sigma", positive = TRUE)
    if (ci == "bootstrap") stop("The bootstrap interval needs `sigma = NULL` (the sample standard deviation).", call. = FALSE)
    if (!is.null(df)) check_number(df, "df", positive = TRUE)
  } else {
    df <- n - 1
  }
  sw <- if (is.null(sigma)) s else sigma

  idx <- function(lsl, usl, sd) {
    cpl <- if (is.null(lsl)) NA_real_ else (m - lsl) / (3 * sd)
    cpu <- if (is.null(usl)) NA_real_ else (usl - m) / (3 * sd)
    cp <- if (is.null(lsl) || is.null(usl)) NA_real_ else (usl - lsl) / (6 * sd)
    c(cp = cp, cpl = cpl, cpu = cpu, cpk = min(cpl, cpu, na.rm = TRUE))
  }
  est <- idx(lsl, usl, sw)
  z <- stats::qnorm(1 - (1 - conf) / 2)
  lower <- upper <- rep(NA_real_, 4); names(lower) <- names(upper) <- names(est)
  lower_one <- NA_real_
  se_k <- function(k) sqrt(1 / (9 * n) + k^2 / (2 * df))
  if (ci == "normal" && !is.null(df)) {
    for (nm in c("cpl", "cpu", "cpk")) {
      if (!is.na(est[[nm]])) {
        lower[[nm]] <- est[[nm]] - z * se_k(est[[nm]]); upper[[nm]] <- est[[nm]] + z * se_k(est[[nm]])
      }
    }
    if (!is.na(est[["cp"]])) {
      lower[["cp"]] <- est[["cp"]] * sqrt(stats::qchisq((1 - conf) / 2, df) / df)
      upper[["cp"]] <- est[["cp"]] * sqrt(stats::qchisq(1 - (1 - conf) / 2, df) / df)
    }
    lower_one <- est[["cpk"]] - stats::qnorm(conf) * se_k(est[["cpk"]])
  } else if (ci == "bootstrap") {
    B <- as.integer(B)
    if (B < 100) stop("`B` must be at least 100.", call. = FALSE)
    xb <- matrix(sample(x, n * B, replace = TRUE), n, B)
    mb <- colMeans(xb); sb <- sqrt((colSums(xb^2) - n * mb^2) / (n - 1))
    bs <- list(
      cpl = if (is.null(lsl)) NULL else (mb - lsl) / (3 * sb),
      cpu = if (is.null(usl)) NULL else (usl - mb) / (3 * sb),
      cp = if (is.null(lsl) || is.null(usl)) NULL else (usl - lsl) / (6 * sb))
    bs$cpk <- do.call(pmin, Filter(Negate(is.null), list(bs$cpl, bs$cpu)))
    a <- (1 - conf) / 2
    for (nm in names(bs)) if (!is.null(bs[[nm]])) {
      q <- stats::quantile(bs[[nm]], c(a, 1 - a), names = FALSE)
      lower[[nm]] <- q[1]; upper[[nm]] <- q[2]
    }
    lower_one <- stats::quantile(bs$cpk, 1 - conf, names = FALSE)
  }
  ppk <- idx(lsl, usl, s)
  cpm <- if (is.null(target) || is.null(lsl) || is.null(usl)) NA_real_ else (usl - lsl) / (6 * sqrt(sw^2 + (m - target)^2))
  indices <- tibble::tibble(
    index = c("Cp", "Cpl", "Cpu", "Cpk", "Pp", "Ppk", "Cpm"),
    estimate = c(est[["cp"]], est[["cpl"]], est[["cpu"]], est[["cpk"]], ppk[["cp"]], ppk[["cpk"]], cpm),
    lower = c(lower[["cp"]], lower[["cpl"]], lower[["cpu"]], lower[["cpk"]], rep(NA_real_, 3)),
    upper = c(upper[["cp"]], upper[["cpl"]], upper[["cpu"]], upper[["cpk"]], rep(NA_real_, 3))
  )
  indices <- indices[!is.na(indices$estimate), ]
  out_obs <- mean((!is.null(lsl) & x < (if (is.null(lsl)) -Inf else lsl)) | (!is.null(usl) & x > (if (is.null(usl)) Inf else usl)))
  exp_out <- (if (is.null(lsl)) 0 else stats::pnorm(lsl, m, sw)) + (if (is.null(usl)) 0 else stats::pnorm(usl, m, sw, lower.tail = FALSE))
  structure(list(indices = indices, n = n, mean = m, sd = s, sigma = sw, lower_one_sided = lower_one,
                 observed_out = out_obs, expected_out = exp_out, dpmo = 1e6 * exp_out, lsl = lsl, usl = usl,
                 conf = conf, ci = ci),
            class = "capability")
}

#' @export
print.capability <- function(x, digits = 3, ...) {
  cat("Process capability: n = ", x$n, ", mean = ", format(round(x$mean, digits)), ", SD = ", format(round(x$sd, digits)),
      if (!isTRUE(all.equal(x$sigma, x$sd))) paste0(", within SD = ", format(round(x$sigma, digits))), "\n",
      "Specification: ", if (!is.null(x$lsl)) paste0("LSL = ", x$lsl) else "no LSL", ", ",
      if (!is.null(x$usl)) paste0("USL = ", x$usl) else "no USL", "\n\n", sep = "")
  d <- as.data.frame(x$indices)
  d[-1] <- lapply(d[-1], function(v) signif(v, digits))
  print(d, row.names = FALSE)
  if (x$ci != "none") {
    cat("\nIntervals: ", format(100 * x$conf), "% ", if (x$ci == "normal") "(Bissell normal approximation; exact for Cp)" else "(percentile bootstrap)",
        "; one-sided lower bound for Cpk = ", format(round(x$lower_one_sided, digits)), "\n", sep = "")
  }
  cat("Outside the specification: observed ", format(round(100 * x$observed_out, 2)), "%, expected under normality ",
      format(round(100 * x$expected_out, 3)), "% (", format(round(x$dpmo, 1)), " per million)\n", sep = "")
  invisible(x)
}

#' Simulate the sampling distribution of the capability indices
#'
#' Draws `nsim` samples of size `n` from a normal process with mean `mu` and standard deviation `sigma`,
#' computes the capability indices \eqn{C_p} and \eqn{C_{pk}} of each sample (with the sample mean and SD), and the
#' Bissell normal-approximation interval and one-sided lower bound for \eqn{C_{pk}} (see [capability()]). The spread of
#' the estimates shows how uncertain a capability index from a small sample is; the interval coverage shows how
#' reliable the approximation is. Use [set.seed()] before the call for a reproducible result.
#'
#' The true indices of the process are in the attributes `true_cp` and `true_cpk`.
#'
#' @param n Sample size (at least 3).
#' @param mu,sigma Mean and standard deviation of the process (`sigma` positive).
#' @param lsl,usl Specification limits.
#' @param nsim Number of simulated samples. Default 10000.
#' @param conf Confidence level of the intervals. Default 0.95.
#' @return A tibble with one row per sample and columns `cp`, `cpk`, `lower`, `upper` (two-sided interval for Cpk), and
#'   `lower_one` (one-sided lower bound), with attributes `true_cp` and `true_cpk`.
#' @seealso [capability()]
#' @examples
#' set.seed(1)
#' # a process with true Cpk = 1.33 (centered), samples of 30
#' s <- sim_capability(30, mu = 500, sigma = 25 / (3 * 1.33), lsl = 475, usl = 525, nsim = 2000)
#' attr(s, "true_cpk")
#' mean(s$cpk >= 1.4)     # how often a supplier at the contract boundary reports 1.4 or more
#' mean(s$lower < 1.33 & s$upper > 1.33)  # coverage of the 95% interval
#' @export
sim_capability <- function(n, mu, sigma, lsl, usl, nsim = 10000, conf = 0.95) {
  check_count(n, "n", minimum = 3); check_count(nsim, "nsim", minimum = 1)
  check_finite(mu, "mu"); check_number(sigma, "sigma", positive = TRUE)
  check_finite(lsl, "lsl"); check_finite(usl, "usl")
  if (lsl >= usl) stop("`lsl` must be below `usl`.", call. = FALSE)
  check_probability(conf, "conf")
  x <- matrix(stats::rnorm(n * nsim, mu, sigma), n, nsim)
  m <- colMeans(x); s <- sqrt((colSums(x^2) - n * m^2) / (n - 1))
  cpk <- pmin(m - lsl, usl - m) / (3 * s)
  se <- sqrt(1 / (9 * n) + cpk^2 / (2 * (n - 1)))
  z <- stats::qnorm(1 - (1 - conf) / 2)
  out <- tibble::tibble(cp = (usl - lsl) / (6 * s), cpk = cpk, lower = cpk - z * se, upper = cpk + z * se,
                        lower_one = cpk - stats::qnorm(conf) * se)
  attr(out, "true_cp") <- (usl - lsl) / (6 * sigma)
  attr(out, "true_cpk") <- min(mu - lsl, usl - mu) / (3 * sigma)
  out
}

#' Defects per million and sigma level
#'
#' `dpmo()` gives the expected defects per million opportunities beyond the **nearer** specification limit of a
#' normal process with capability index `cpk`, allowing a shift of the process mean toward that limit (in units
#' of the standard deviation): \eqn{10^6\,[1 - \Phi(3C_{pk} - \text{shift})]}{1e6 (1 - Phi(3 Cpk - shift))}. `sigma_level()` converts defects per million back to a sigma level,
#' \eqn{\Phi^{-1}(1 - \text{DPMO}/10^6) + \text{shift}}{qnorm(1 - DPMO/1e6) + shift}, with the conventional shift of 1.5.
#'
#' With a shift of 1.5, a process with \eqn{C_{pk} = 2} (the "six sigma" process) has 3.4 defects per million.
#'
#' @param cpk Capability index (a number or a vector).
#' @param shift Shift of the mean toward the nearer limit, in standard deviations. Default 0.
#' @param dpmo Defects per million (a number or a vector between 0 and 1e6, exclusive).
#' @return A numeric vector.
#' @examples
#' dpmo(c(1, 1.33, 1.67, 2))              # no shift
#' dpmo(2, shift = 1.5)                     # 3.4: the six sigma convention
#' sigma_level(c(66807, 3.4))               # 3 and 6
#' @export
dpmo <- function(cpk, shift = 0) {
  if (!is.numeric(cpk) || anyNA(cpk)) stop("`cpk` must be numeric with no missing values.", call. = FALSE)
  check_finite(shift, "shift")
  1e6 * stats::pnorm(3 * cpk - shift, lower.tail = FALSE)
}

#' @rdname dpmo
#' @export
sigma_level <- function(dpmo, shift = 1.5) {
  if (!is.numeric(dpmo) || anyNA(dpmo) || any(dpmo <= 0 | dpmo >= 1e6)) stop("`dpmo` must be between 0 and 1,000,000 (exclusive).", call. = FALSE)
  check_finite(shift, "shift")
  stats::qnorm(1 - dpmo / 1e6) + shift
}

check_finite <- function(x, name) {
  if (!(is.numeric(x) && length(x) == 1 && is.finite(x))) stop("`", name, "` must be a single finite number.", call. = FALSE)
  invisible(x)
}

#' Acceptance of a supplier by a capability index: the OC curve and a sampling plan
#'
#' A supplier is accepted when the capability index \eqn{\hat C_{pk}}{Cpk-hat} of a sample of `n` units is at least
#' `k`. `cpk_oc()` gives the probability of acceptance as a function of the true \eqn{C_{pk}}{Cpk}: the **OC curve** of
#' the rule (as for the sampling plans of Topic 3), using the large-sample normal approximation of Bissell,
#' \eqn{\hat C_{pk} \approx N(C_{pk},\, 1/(9n) + C_{pk}^2/(2(n-1)))}{Cpk-hat ~ N(Cpk, 1/(9n) + Cpk^2/(2(n - 1)))}. `capability_plan()`
#' finds the smallest `n`, and the cut-off `k`, such that a supplier with true \eqn{C_{pk}}{Cpk} equal to `cpk_good` is
#' rejected with probability `alpha` (the producer's risk) and one with `cpk_bad` is accepted with probability at most
#' `beta` (the consumer's risk).
#'
#' The approximation is rough for small samples: at `n = 30` an acceptance probability can differ from the simulated one by
#' about 0.06 (the sampling distribution of the sample Cpk is skewed). Check a plan with [sim_capability()].
#'
#' @param cpk True capability index (a number or a vector).
#' @param n Sample size.
#' @param k Cut-off: accept if the sample Cpk is at least `k`.
#' @param cpk_good,cpk_bad The good and the bad true Cpk, with `cpk_bad < cpk_good`.
#' @param alpha,beta Producer's risk and consumer's risk.
#' @param max_n Largest sample size searched. Default 5000.
#' @return `cpk_oc()` returns a numeric vector of acceptance probabilities. `capability_plan()` returns a one-row tibble
#'   with `n`, `k`, `producer_risk` and `consumer_risk` (both at the approximation).
#' @seealso [capability()], [sim_capability()]
#' @examples
#' cpk_oc(c(1, 1.33, 1.67), n = 30, k = 1.4)
#' capability_plan(cpk_good = 1.67, cpk_bad = 1.33)
#' @export
cpk_oc <- function(cpk, n, k) {
  if (!is.numeric(cpk) || anyNA(cpk) || any(cpk <= 0)) stop("`cpk` must be positive numbers.", call. = FALSE)
  check_count(n, "n", minimum = 3); check_finite(k, "k")
  se <- sqrt(1 / (9 * n) + cpk^2 / (2 * (n - 1)))
  stats::pnorm((cpk - k) / se)
}

#' @rdname cpk_oc
#' @export
capability_plan <- function(cpk_good, cpk_bad, alpha = 0.05, beta = 0.10, max_n = 5000) {
  check_number(cpk_good, "cpk_good", positive = TRUE); check_number(cpk_bad, "cpk_bad", positive = TRUE)
  if (cpk_bad >= cpk_good) stop("`cpk_bad` must be below `cpk_good`.", call. = FALSE)
  check_probability(alpha, "alpha"); check_probability(beta, "beta"); check_count(max_n, "max_n", minimum = 4)
  for (n in 4:max_n) {
    se_good <- sqrt(1 / (9 * n) + cpk_good^2 / (2 * (n - 1)))
    k <- cpk_good - stats::qnorm(1 - alpha) * se_good
    if (cpk_oc(cpk_bad, n, k) <= beta) {
      return(tibble::tibble(n = n, k = k, producer_risk = 1 - cpk_oc(cpk_good, n, k), consumer_risk = cpk_oc(cpk_bad, n, k)))
    }
  }
  stop("No plan with at most `max_n` units meets the risks.", call. = FALSE)
}
