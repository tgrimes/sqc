#' Variance components of a one-way random effects model (ANOVA estimators)
#'
#' Fits the one-way random effects model
#' \eqn{y_{ij} = \mu + a_i + e_{ij}}{y_ij = mu + a_i + e_ij}, where the group effects
#' \eqn{a_i \sim N(0, \sigma_a^2)}{a_i ~ N(0, sigma_a^2)} and the errors
#' \eqn{e_{ij} \sim N(0, \sigma_e^2)}{e_ij ~ N(0, sigma_e^2)} are independent, and estimates the two
#' variance components by the **ANOVA (method of moments) estimators**
#' \deqn{\hat\sigma_e^2 = MSW, \qquad \hat\sigma_a^2 = (MSB - MSW)/n_0,}{sigma_e^2-hat = MSW, sigma_a^2-hat = (MSB - MSW)/n0,}
#' where \eqn{MSB} and \eqn{MSW} are the between-group and within-group mean squares and
#' \eqn{n_0 = (N - \sum n_i^2/N)/(a - 1)}{n0 = (N - sum(n_i^2)/N)/(a - 1)} equals the common group size
#' \eqn{n} when the design is balanced.
#'
#' The **intraclass correlation** is \eqn{\text{ICC} = \sigma_a^2/(\sigma_a^2 + \sigma_e^2)}{ICC = sigma_a^2/(sigma_a^2 + sigma_e^2)}:
#' the share of the total variance that is due to the groups, and also the correlation between two
#' observations from the same group. It is estimated by \eqn{(F - 1)/(F + n_0 - 1)}{(F - 1)/(F + n0 - 1)}
#' with \eqn{F = MSB/MSW}.
#'
#' The ANOVA estimator of \eqn{\sigma_a^2}{sigma_a^2} is unbiased but can be **negative**, when
#' \eqn{MSB < MSW}. The estimate is then reported as it is in the `variance` column of
#' `components` (a warning flag is in `negative`), and the standard deviation, the percent of the
#' total, and the ICC use 0 instead. A negative estimate says that the data cannot separate the
#' group component from zero, not that the variance is negative.
#'
#' The confidence interval for the ICC is the exact F-based interval for a balanced design and an
#' approximation (with \eqn{n_0}{n0}) for an unbalanced one. The interval for \eqn{\sigma_e}{sigma_e} is the
#' exact chi-square interval. Limits of the ICC interval are truncated to 0 and 1.
#'
#' No current course topic calls `vc_oneway()`; it is kept for reference (Topic 11's technical details on random
#' blocks estimate the same intraclass correlation from the block ANOVA directly).
#'
#' @param data A matrix with one row per group and one column per replicate (a balanced design), or a
#'   data frame with one row per observation (use `value` and `group`). Groups may have different
#'   sizes in a data frame; each group needs at least 1 observation, at least 2 groups are needed,
#'   and at least one group must have 2 or more.
#' @param value,group For a data frame: the column holding the measurements and the column
#'   identifying the group (bare names, as in dplyr).
#' @param conf Confidence level of the intervals. Default 0.95.
#' @return An object of class `vc_oneway`, a list with elements `anova` (a tibble with the
#'   source, degrees of freedom, sum of squares, mean square, F, and p-value), `components` (a
#'   tibble with the variance, standard deviation, and percent of the total variance of the
#'   between-group and within-group components), `icc` (estimate and confidence limits), `sd_within_ci`
#'   (confidence interval for \eqn{\sigma_e}{sigma_e}), `negative`, `n0`, `a` (number of groups), `N`
#'   (number of observations), `balanced`, and `conf`.
#' @seealso [sim_vc()] to simulate the sampling distribution of the estimators, and [gage_rr()] for the crossed
#'   two-factor model.
#' @examples
#' # The assay_study precision study: ten samples measured six times each
#' # (3 technicians x 2 replicates)
#' vc_oneway(assay_study, value = hba1c, group = sample)
#'
#' # A matrix with one row per group
#' set.seed(1)
#' y <- matrix(rnorm(8 * 4, mean = 100, sd = 2), ncol = 4) + rnorm(8, sd = 3)
#' vc_oneway(y)
#' @export
vc_oneway <- function(data, value = NULL, group = NULL, conf = 0.95) {
  value <- rlang::enquo(value)
  group <- rlang::enquo(group)
  if (!(is.numeric(conf) && length(conf) == 1 && is.finite(conf) && conf > 0 && conf < 1)) {
    stop("`conf` must be a single number between 0 and 1.", call. = FALSE)
  }
  vg <- value_group(data, value, group)
  y <- vg$y
  g <- vg$g
  a <- nlevels(g)
  N <- length(y)
  if (a < 2) stop("At least 2 groups are needed.", call. = FALSE)
  sizes <- as.vector(table(g))
  if (N - a < 1) stop("At least one group needs 2 or more observations.", call. = FALSE)

  group_mean <- tapply(y, g, mean)
  grand_mean <- mean(y)
  ssb <- sum(sizes * (group_mean - grand_mean)^2)
  ssw <- sum((y - group_mean[g])^2)
  df_b <- a - 1
  df_w <- N - a
  msb <- ssb / df_b
  msw <- ssw / df_w
  f <- msb / msw
  p_value <- stats::pf(f, df_b, df_w, lower.tail = FALSE)
  n0 <- (N - sum(sizes^2) / N) / (a - 1)
  balanced <- length(unique(sizes)) == 1

  s2_e <- msw
  s2_a <- (msb - msw) / n0
  negative <- s2_a < 0
  var_pos <- c(max(0, s2_a), s2_e)
  total <- sum(var_pos)
  components <- tibble::tibble(
    component = c("between groups", "within groups", "total"),
    variance = c(s2_a, s2_e, total),
    sd = sqrt(pmax(0, c(s2_a, s2_e, total))),
    percent = 100 * c(var_pos, total) / total
  )

  alpha <- 1 - conf
  icc_of <- function(ff) pmin(1, pmax(0, (ff - 1) / (ff + n0 - 1)))
  f_lower <- f / stats::qf(1 - alpha / 2, df_b, df_w)
  f_upper <- f * stats::qf(1 - alpha / 2, df_w, df_b)
  icc <- c(estimate = icc_of(f), lower = icc_of(f_lower), upper = icc_of(f_upper))
  sd_within_ci <- c(estimate = sqrt(s2_e),
                    lower = sqrt(ssw / stats::qchisq(1 - alpha / 2, df_w)),
                    upper = sqrt(ssw / stats::qchisq(alpha / 2, df_w)))

  structure(
    list(
      anova = tibble::tibble(
        source = c("Between groups", "Within groups"),
        df = c(df_b, df_w), ss = c(ssb, ssw), ms = c(msb, msw),
        f = c(f, NA_real_), p = c(p_value, NA_real_)
      ),
      components = components, icc = icc, sd_within_ci = sd_within_ci,
      negative = negative, n0 = n0, a = a, N = N, balanced = balanced, conf = conf
    ),
    class = "vc_oneway"
  )
}

# Measurements and group labels from a matrix (one row per group) or a data frame.
value_group <- function(data, value, group) {
  if (is.data.frame(data)) {
    if (rlang::quo_is_null(value) || rlang::quo_is_null(group)) {
      stop("For a data frame, supply both `value` and `group` (the measurement and the group columns). ",
           "If each row is already a group (wide format), convert it with `as.matrix()`.", call. = FALSE)
    }
    y <- rlang::eval_tidy(value, data)
    g <- rlang::eval_tidy(group, data)
    g <- factor(g, levels = unique(g))
  } else if (is.matrix(data)) {
    if (!rlang::quo_is_null(value) || !rlang::quo_is_null(group)) {
      stop("`value` and `group` are only used with a data frame.", call. = FALSE)
    }
    y <- as.vector(t(data))
    g <- factor(rep(seq_len(nrow(data)), each = ncol(data)))
  } else {
    stop("`data` must be a matrix (one row per group) or a data frame.", call. = FALSE)
  }
  if (!is.numeric(y) || anyNA(y) || anyNA(g)) {
    stop("Measurements must be numeric, and measurements and groups must have no missing values.", call. = FALSE)
  }
  list(y = y, g = g)
}

#' @export
print.vc_oneway <- function(x, digits = 4, ...) {
  cat("One-way random effects model: ", x$a, " groups, ", x$N, " observations",
      if (!x$balanced) " (unbalanced)", "\n\n", sep = "")
  an <- as.data.frame(x$anova)
  an$source <- as.character(an$source)
  print(format_table(an, digits), row.names = FALSE)
  cat("\n")
  comp <- as.data.frame(x$components)
  comp$percent <- round(comp$percent, 1)
  print(format_table(comp, digits), row.names = FALSE)
  if (x$negative) {
    cat("\nThe ANOVA estimate of the between-group variance is negative; the sd, percent, and ICC use 0.\n")
  }
  cat("\nICC = ", format(round(x$icc[["estimate"]], digits), nsmall = 0),
      "  (", format(100 * x$conf), "% CI ", format(round(x$icc[["lower"]], digits)), " to ",
      format(round(x$icc[["upper"]], digits)), ")\n", sep = "")
  cat("SD within groups = ", format(round(x$sd_within_ci[["estimate"]], digits)),
      "  (", format(100 * x$conf), "% CI ", format(round(x$sd_within_ci[["lower"]], digits)), " to ",
      format(round(x$sd_within_ci[["upper"]], digits)), ")\n", sep = "")
  invisible(x)
}

format_table <- function(d, digits) {
  d[] <- lapply(d, function(col) if (is.numeric(col)) signif(col, digits) else col)
  d
}

#' Simulate the ANOVA estimators of a one-way random effects model
#'
#' Generates `nsim` balanced data sets from the one-way random effects model
#' \eqn{y_{ij} = a_i + e_{ij}}{y_ij = a_i + e_ij} with \eqn{a_i \sim N(0, \sigma_a^2)}{a_i ~ N(0, sigma_a^2)},
#' \eqn{e_{ij} \sim N(0, \sigma_e^2)}{e_ij ~ N(0, sigma_e^2)}, `a` groups, and `n` observations per
#' group, and computes the ANOVA estimators for each data set (see [vc_oneway()]). The averages
#' of the estimates of \eqn{\sigma_a^2}{sigma_a^2} and \eqn{\sigma_e^2}{sigma_e^2} are close to the truth (the
#' estimators are unbiased), but their spread shows how precisely a study of this size estimates
#' them, and the share of negative estimates of \eqn{\sigma_a^2}{sigma_a^2} shows how often the data cannot
#' separate a small group component from zero.
#'
#' Use [set.seed()] before the call for a reproducible result.
#'
#' @param a Number of groups (at least 2).
#' @param n Number of observations per group (at least 2).
#' @param sigma_a,sigma_e Standard deviations of the group effects (non-negative) and of the errors
#'   (positive).
#' @param nsim Number of simulated data sets. Default 10000.
#' @return A tibble with one row per data set and columns `msb`, `msw` (the mean squares), `f`
#'   (their ratio), `sigma2_a` (the ANOVA estimate of \eqn{\sigma_a^2}{sigma_a^2}, which can be negative),
#'   `sigma2_e`, and `icc` (the ANOVA estimate of the ICC, which is negative when `f < 1`).
#' @seealso [vc_oneway()]
#' @examples
#' set.seed(1)
#' s <- sim_vc(a = 10, n = 6, sigma_a = 1.2, sigma_e = 0.25, nsim = 2000)
#' colMeans(s[c("sigma2_a", "sigma2_e", "icc")])   # truth: 1.44, 0.0625, 0.96
#' mean(s$sigma2_a < 0)                              # share of negative estimates
#' @export
sim_vc <- function(a, n, sigma_a, sigma_e, nsim = 10000) {
  check_count(a, "a", minimum = 2)
  check_count(n, "n", minimum = 2)
  check_count(nsim, "nsim", minimum = 1)
  if (!(is.numeric(sigma_a) && length(sigma_a) == 1 && is.finite(sigma_a) && sigma_a >= 0)) {
    stop("`sigma_a` must be a single non-negative finite number.", call. = FALSE)
  }
  check_number(sigma_e, "sigma_e", positive = TRUE)

  group_effect <- matrix(stats::rnorm(a * nsim, sd = sigma_a), a, nsim)
  y <- matrix(stats::rnorm(n * a * nsim, sd = sigma_e), n, a * nsim) + rep(as.vector(group_effect), each = n)
  group_mean <- colMeans(y)                                    # one mean per group per data set
  ssw <- colSums(matrix(colSums((y - rep(group_mean, each = n))^2), a, nsim))
  gm <- matrix(group_mean, a, nsim)
  ssb <- n * colSums((gm - rep(colMeans(gm), each = a))^2)
  msb <- ssb / (a - 1)
  msw <- ssw / (a * (n - 1))
  tibble::tibble(msb = msb, msw = msw, f = msb / msw, sigma2_a = (msb - msw) / n,
                 sigma2_e = msw, icc = (msb / msw - 1) / (msb / msw + n - 1))
}

check_count <- function(x, name, minimum = 1) {
  if (!(is.numeric(x) && length(x) == 1 && is.finite(x) && x == round(x) && x >= minimum)) {
    stop("`", name, "` must be a single whole number of at least ", minimum, ".", call. = FALSE)
  }
  invisible(x)
}
