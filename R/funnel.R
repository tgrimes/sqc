#' Funnel plot for comparing institutions
#'
#' A **funnel plot** compares institutions that differ in size without ranking them. Each institution has a count \eqn{x_i} and a
#' size \eqn{s_i}, and its plotted value \eqn{x_i/s_i} is compared with a **target** \eqn{\theta_0}. The variability of the
#' value shrinks as the size grows, so the control limits form a funnel around the target: the limits at the 95% and 99.8%
#' levels (about 2 and 3 standard deviations) are the Poisson or binomial limits for an institution of that size.
#'
#' * `type = "ratio"`: `count` is the observed number of events and `size` the **expected** number (the standardized
#'   infection ratio is `count / size`). Under the target ratio \eqn{\theta_0 = 1} the count is Poisson with mean \eqn{\theta_0 s_i}.
#' * `type = "proportion"`: `count` is the number of events and `size` the number of cases. The target is the proportion
#'   \eqn{\theta_0} (default: the pooled proportion) and the count is binomial\eqn{(s_i, \theta_0)}.
#'
#' Because the counts are discrete, institutions are judged with the **mid-p** z-score
#' \eqn{z_i = \Phi^{-1}\{P(X < x_i) + \tfrac12 P(X = x_i)\}}, and an institution is outside the limits at level \eqn{L} when
#' \eqn{|z_i| > \Phi^{-1}\{1 - (1 - L)/2\}} (1.96 for 95% and 3.09 for 99.8%). With a large number of institutions some will be
#' outside the 95% limits by chance (5% of them) and very few outside the 99.8% limits (0.2%).
#'
#' **Overdispersion.** If the institutions differ by more than chance (for example because of unmeasured case mix), too many
#' are flagged. Following Spiegelhalter (2005), a dispersion factor \eqn{\phi} is estimated from the **winsorized** z-scores (the lower
#' and upper `winsor` fraction replaced by the quantiles; the mean of the squared z-scores is divided by its value for normal data, so
#' that \eqn{\phi} is about 1 without overdispersion). `overdispersion = "multiplicative"` widens the limits by \eqn{\sqrt\phi} on the
#' z-scale. `overdispersion = "additive"` estimates a between-institution variance \eqn{\tau^2} by the DerSimonian-Laird method and,
#' for ratios, uses the negative binomial distribution (Poisson counts whose true ratios vary with variance \eqn{\tau^2}) for the
#' z-scores and limits; for proportions it uses the normal limits \eqn{\theta_0 \pm z_L\sqrt{\theta_0(1 - \theta_0)/s_i + \tau^2}}. Either adjustment is only
#' applied if \eqn{\phi > 1}.
#'
#' @param count Non-negative whole numbers: the observed counts.
#' @param size Positive numbers: the expected counts (`type = "ratio"`) or the numbers of cases (`type = "proportion"`).
#' @param type `"ratio"` (default) or `"proportion"`.
#' @param target The target ratio or proportion. Default 1 for `"ratio"` and the pooled proportion for `"proportion"`.
#' @param levels The two confidence levels of the limits. Default `c(0.95, 0.998)`.
#' @param overdispersion `"none"` (default), `"multiplicative"`, or `"additive"`.
#' @param winsor Fraction winsorized at each end when estimating \eqn{\phi}. Default 0.1.
#' @param label Optional labels (for example hospital names), used by [autoplot.sqc_funnel()].
#' @return An object of class `sqc_funnel`: a list with `statistics` (a tibble with `label`, `count`, `size`, `value`, `z`, the
#'   adjusted `z_adj`, the columns `outside_95` and `outside_998` coded `""`, `"high"`, or `"low"` for the two levels, and
#'   `status`), `curves` (the limits on a grid of sizes), `target`, `type`, `levels`, `overdispersion`, `phi`, and `tau2`.
#' @references Spiegelhalter, D. J. (2005). Funnel plots for comparing institutional performance. *Statistics in Medicine*, 24, 1185-1202.
#' @seealso [sim_funnel()], [count_chart()]
#' @examples
#' f <- funnel_plot(hospital_infections$observed, hospital_infections$predicted,
#'                  label = hospital_infections$hospital)
#' f
#' autoplot(f)
#' @export
funnel_plot <- function(count, size, type = c("ratio", "proportion"), target = NULL, levels = c(0.95, 0.998),
                        overdispersion = c("none", "multiplicative", "additive"), winsor = 0.1, label = NULL) {
  type <- match.arg(type); overdispersion <- match.arg(overdispersion)
  fit <- funnel_fit(count, size, type, target, levels, overdispersion, winsor)
  m <- length(count)
  if (is.null(label)) label <- as.character(seq_len(m))
  if (length(label) != m) stop("`label` must have one value per institution.", call. = FALSE)
  side <- function(j) ifelse(fit$z_adj > fit$zcrit[j], "high", ifelse(fit$z_adj < -fit$zcrit[j], "low", ""))
  o1 <- side(1); o2 <- side(2)
  status <- ifelse(o2 != "", paste0(levels[2] * 100, "% limits"), ifelse(o1 != "", paste0(levels[1] * 100, "% limits"), "inside"))
  stat <- tibble::tibble(label = label, count = count, size = size, value = count / size, z = fit$z, z_adj = fit$z_adj,
                         outside_95 = o1, outside_998 = o2, status = status)
  curves <- funnel_curves(range(size), type, fit, levels)
  structure(list(statistics = stat, curves = curves, target = fit$target, type = type, levels = levels, overdispersion = overdispersion,
                 phi = fit$phi, tau2 = fit$tau2), class = "sqc_funnel")
}

# z-scores, overdispersion estimates, and the critical values.
funnel_fit <- function(count, size, type, target, levels, overdispersion, winsor) {
  if (!is.numeric(count) || anyNA(count) || any(count < 0) || any(count != round(count)) || length(count) < 3) {
    stop("`count` must be at least 3 non-negative whole numbers with no missing values.", call. = FALSE)
  }
  if (!is.numeric(size) || anyNA(size) || any(size <= 0) || length(size) != length(count)) {
    stop("`size` must be positive numbers, one for each count.", call. = FALSE)
  }
  if (type == "proportion" && (any(count > size) || any(size != round(size)))) {
    stop("For proportions, `size` must be whole numbers of cases and `count` cannot exceed it.", call. = FALSE)
  }
  if (!(is.numeric(levels) && length(levels) == 2 && all(levels > 0 & levels < 1) && levels[1] < levels[2])) {
    stop("`levels` must be two probabilities, the smaller first.", call. = FALSE)
  }
  if (!(is.numeric(winsor) && length(winsor) == 1 && winsor >= 0 && winsor < 0.5)) stop("`winsor` must be between 0 and 0.5.", call. = FALSE)
  if (is.null(target)) target <- if (type == "ratio") 1 else sum(count) / sum(size)
  check_number(target, "target", positive = TRUE)
  if (type == "proportion" && target >= 1) stop("`target` must be a proportion below 1.", call. = FALSE)
  zcrit <- stats::qnorm(1 - (1 - levels) / 2)
  z <- zmid(count, size, target, type)
  zw <- pmin(pmax(z, stats::quantile(z, winsor, names = FALSE)), stats::quantile(z, 1 - winsor, names = FALSE))
  # The winsorized mean of z^2 is below 1 even without overdispersion; divide by its value for N(0, 1) data so that phi is about 1 then.
  cw <- if (winsor > 0) { q <- stats::qnorm(1 - winsor); 2 * ((stats::pnorm(q) - 0.5) - q * stats::dnorm(q) + q^2 * (1 - stats::pnorm(q))) } else 1
  phi <- mean(zw^2) / cw; tau2 <- 0
  z_adj <- z; phi_used <- 1
  if (overdispersion == "multiplicative" && phi > 1) {
    phi_used <- phi; z_adj <- z / sqrt(phi)
  } else if (overdispersion == "additive" && phi > 1) {
    s2 <- if (type == "ratio") target / size else target * (1 - target) / size
    w <- 1 / s2
    tau2 <- max(0, (length(count) * phi - (length(count) - 1)) / (sum(w) - sum(w^2) / sum(w)))
    z_adj <- if (tau2 == 0) z else if (type == "ratio") zmid(count, size, target, type, tau2) else (count / size - target) / sqrt(s2 + tau2)
  }
  list(z = z, z_adj = z_adj, zcrit = zcrit, phi = phi, phi_used = phi_used, tau2 = tau2, target = target, overdispersion = overdispersion)
}

# Mid-p z-score of a count under the Poisson (ratio) or binomial (proportion) null.
zmid <- function(count, size, target, type, tau2 = 0) {
  if (type == "ratio" && tau2 > 0) {
    # additive between-institution variation: the counts are negative binomial with variance mu + (tau2 / target^2) mu^2
    k <- target^2 / tau2; mu <- target * size
    below <- stats::pnbinom(count - 1, size = k, mu = mu); at <- stats::dnbinom(count, size = k, mu = mu)
  } else if (type == "ratio") {
    below <- stats::ppois(count - 1, target * size); at <- stats::dpois(count, target * size)
  } else {
    below <- stats::pbinom(count - 1, size, target); at <- stats::dbinom(count, size, target)
  }
  stats::qnorm(pmin(pmax(below + 0.5 * at, 1e-300), 1 - 1e-16))
}

# The limits (as plotted values) on a grid of sizes, for both levels.
funnel_curves <- function(rng, type, fit, levels) {
  grid <- if (type == "ratio") seq(rng[1] * 0.9, rng[2] * 1.05, length.out = 200) else unique(round(seq(rng[1] * 0.9, rng[2] * 1.05, length.out = 200)))
  grid <- pmax(grid, if (type == "ratio") 0.05 else 1)
  one <- function(level_index) {
    zc <- fit$zcrit[level_index]
    if (fit$overdispersion == "additive" && fit$tau2 > 0 && type == "proportion") {
      s2 <- fit$target * (1 - fit$target) / grid
      w <- zc * sqrt(s2 + fit$tau2)
      return(tibble::tibble(size = grid, lower = pmax(0, fit$target - w), upper = fit$target + w))
    }
    crit <- zc * sqrt(fit$phi_used)
    t2 <- if (fit$overdispersion == "additive" && type == "ratio") fit$tau2 else 0
    up <- lo <- numeric(length(grid))
    for (i in seq_along(grid)) {
      s <- grid[i]
      qf <- function(pr) {
        if (type == "ratio" && t2 > 0) stats::qnbinom(pr, size = fit$target^2 / t2, mu = fit$target * s)
        else if (type == "ratio") stats::qpois(pr, fit$target * s) else stats::qbinom(pr, s, fit$target)
      }
      q <- qf(stats::pnorm(crit))
      xs <- max(0, q - 3):(q + 6)
      zz <- zmid(xs, s, fit$target, type, t2)
      up[i] <- xs[which(zz > crit)[1]] / s
      q2 <- qf(stats::pnorm(-crit))
      xl <- max(0, q2 - 6):(q2 + 3)
      zl <- zmid(xl, s, fit$target, type, t2)
      lo[i] <- if (any(zl < -crit)) xl[max(which(zl < -crit))] / s else NA_real_
    }
    tibble::tibble(size = grid, lower = lo, upper = up)
  }
  lab <- paste0(levels * 100, "%")
  a <- one(1); a$level <- lab[1]; b <- one(2); b$level <- lab[2]
  tibble::as_tibble(rbind(a, b))
}

#' @export
print.sqc_funnel <- function(x, ...) {
  s <- x$statistics
  cat("<sqc_funnel> ", nrow(s), " institutions, ", x$type, " funnel plot, target ", format(signif(x$target, 4)), "\n", sep = "")
  cat("  outside the ", x$levels[2] * 100, "% limits: ", sum(s$outside_998 != ""), " (high ", sum(s$outside_998 == "high"), ", low ",
      sum(s$outside_998 == "low"), "); outside the ", x$levels[1] * 100, "% limits: ", sum(s$outside_95 != ""), " (expected by chance about ",
      format(round(nrow(s) * (1 - x$levels[1]), 1)), ")\n", sep = "")
  if (x$overdispersion != "none") {
    cat("  overdispersion (", x$overdispersion, "): phi = ", format(round(x$phi, 2)), if (x$tau2 > 0) paste0(", tau^2 = ", format(signif(x$tau2, 3))), "\n", sep = "")
  } else cat("  dispersion factor phi = ", format(round(x$phi, 2)), " (1 means no excess variation)\n", sep = "")
  invisible(x)
}

#' @export
as.data.frame.sqc_funnel <- function(x, ...) as.data.frame(x$statistics)

#' Draw a funnel plot with ggplot2
#'
#' Plots each institution's value against its size, with the target and the limits at the two levels (dashed for the lower level, solid
#' for the higher). Institutions outside the limits are drawn in color and, if `label_outside` is `TRUE`, labeled.
#'
#' @param object A [funnel_plot()] object.
#' @param title,subtitle Plot title and subtitle.
#' @param xlab,ylab Axis labels.
#' @param label_outside Label the institutions outside the higher-level limits? Default `TRUE`.
#' @param log_x Use a log scale for the size? Default `FALSE`.
#' @param ylim Limits of the vertical axis. Default: from the data, with a margin, so that the very wide limits of the smallest institutions do not
#'   flatten the plot (the funnel is cut at the top).
#' @param ... Not used.
#' @return A ggplot object.
#' @export
autoplot.sqc_funnel <- function(object, title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL, label_outside = TRUE, log_x = FALSE, ylim = NULL, ...) {
  s <- object$statistics; cv <- object$curves
  if (is.null(xlab)) xlab <- if (object$type == "ratio") "Expected count" else "Number of cases"
  if (is.null(ylab)) ylab <- if (object$type == "ratio") "Observed / expected" else "Proportion"
  if (is.null(title)) title <- "Funnel plot"
  lv <- unique(cv$level)
  s$status <- factor(s$status, levels = c("inside", paste0(object$levels[1] * 100, "% limits"), paste0(object$levels[2] * 100, "% limits")))
  g <- ggplot2::ggplot() +
    ggplot2::geom_line(data = cv, ggplot2::aes(x = .data$size, y = .data$upper, linetype = .data$level), na.rm = TRUE) +
    ggplot2::geom_line(data = cv, ggplot2::aes(x = .data$size, y = .data$lower, linetype = .data$level), na.rm = TRUE) +
    ggplot2::geom_hline(yintercept = object$target, color = "grey30") +
    ggplot2::geom_point(data = s, ggplot2::aes(x = .data$size, y = .data$value, color = .data$status), size = 1.9) +
    ggplot2::scale_linetype_manual(values = stats::setNames(c("dashed", "solid"), lv), name = "Limits") +
    ggplot2::scale_color_manual(values = c("grey40", "#f58518", "#e45756"), drop = FALSE, name = NULL) +
    ggplot2::labs(x = xlab, y = ylab, title = title, subtitle = subtitle) +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
  if (label_outside && any(s$outside_998 != "")) {
    g <- g + ggplot2::geom_text(data = s[s$outside_998 != "", ], ggplot2::aes(x = .data$size, y = .data$value, label = .data$label),
                                hjust = -0.2, vjust = 0.5, size = 3.2)
  }
  if (log_x) g <- g + ggplot2::scale_x_log10()
  if (is.null(ylim)) {
    r <- range(s$value); pad <- 0.15 * diff(r)
    ylim <- if (object$type == "ratio") c(0, max(r[2] * 1.1, object$target * 1.8)) else c(max(0, r[1] - pad), min(1, r[2] + pad))
  }
  g + ggplot2::coord_cartesian(ylim = ylim)
}

#' Simulate how often a funnel plot flags institutions
#'
#' Simulates institutions of the given sizes whose true values differ by a **between-institution standard deviation** `tau`
#' (on the log scale for ratios and the logit scale for proportions; `tau = 0` means that all institutions are exactly at the
#' target), builds the funnel plot of each simulated set of counts with [funnel_plot()], and records how many institutions are
#' outside the limits. With `tau = 0` the proportions flagged are the design values (about 5% beyond the 95% limits and 0.2% beyond
#' the 99.8% limits); with `tau > 0` and no adjustment the funnel plot **over-flags**.
#'
#' @inheritParams funnel_plot
#' @param reps Number of simulated sets of institutions. Default 1000.
#' @param tau Between-institution standard deviation. Default 0.
#' @return A tibble with one row per replicate and the numbers of institutions beyond the lower-level limits (`out_95`) and beyond the
#'   higher-level limits (`out_998`), each split into `high` and `low`, and the total number of institutions `m`.
#' @seealso [funnel_plot()]
#' @examples
#' set.seed(1)
#' s <- sim_funnel(hospital_infections$predicted, tau = 0.25, reps = 200)
#' mean(s$out_998) / s$m[1]     # share beyond the 99.8% limits, far above 0.2%
#' @export
sim_funnel <- function(size, type = c("ratio", "proportion"), target = NULL, tau = 0, reps = 1000, levels = c(0.95, 0.998),
                       overdispersion = c("none", "multiplicative", "additive"), winsor = 0.1) {
  type <- match.arg(type); overdispersion <- match.arg(overdispersion)
  if (!is.numeric(size) || anyNA(size) || any(size <= 0) || length(size) < 3) stop("`size` must be at least 3 positive numbers.", call. = FALSE)
  if (!(is.numeric(tau) && length(tau) == 1 && tau >= 0)) stop("`tau` must be a single number of at least 0.", call. = FALSE)
  check_count(reps, "reps", minimum = 1)
  m <- length(size)
  if (type == "proportion" && is.null(target)) stop("`target` (the true proportion) is required for `type = \"proportion\"`.", call. = FALSE)
  th <- if (is.null(target)) 1 else target
  out <- matrix(NA_real_, reps, 4, dimnames = list(NULL, c("out_95_high", "out_95_low", "out_998_high", "out_998_low")))
  for (r in seq_len(reps)) {
    eps <- stats::rnorm(m)
    if (type == "ratio") {
      mu <- th * exp(tau * eps - tau^2 / 2) * size
      cnt <- stats::rpois(m, mu)
    } else {
      p <- stats::plogis(stats::qlogis(th) + tau * eps)
      cnt <- stats::rbinom(m, size, p)
    }
    f <- funnel_fit(cnt, size, type, target, levels, overdispersion, winsor)
    out[r, ] <- c(sum(f$z_adj > f$zcrit[1]), sum(f$z_adj < -f$zcrit[1]), sum(f$z_adj > f$zcrit[2]), sum(f$z_adj < -f$zcrit[2]))
  }
  tibble::tibble(rep = seq_len(reps), m = m, out_95_high = out[, 1], out_95_low = out[, 2], out_95 = out[, 1] + out[, 2],
                 out_998_high = out[, 3], out_998_low = out[, 4], out_998 = out[, 3] + out[, 4])
}
