#' Tabular CUSUM chart
#'
#' Computes the two-sided **tabular CUSUM** of a series of results in time order. With standardized results
#' \eqn{z_t = (x_t - \mu_0)/\sigma} (the target \eqn{\mu_0} and the standard deviation \eqn{\sigma} of the plotted
#' statistic), the upper and lower cumulative sums are
#' \deqn{C_t^+ = \max(0,\; C_{t-1}^+ + z_t - k), \qquad C_t^- = \max(0,\; C_{t-1}^- - z_t - k),}{C+_t = max(0, C+_(t-1) + z_t - k), C-_t = max(0, C-_(t-1) - z_t - k),}
#' starting at the `head_start`, and the chart **signals** when \eqn{C_t^+ > h} (an upward shift) or \eqn{C_t^- > h} (a downward
#' shift). The reference value `k` is half of the shift, in standard deviations, that the chart is designed to detect
#' (\eqn{k = \delta/2}{k = delta/2}), and the decision interval `h` is chosen to give the in-control average run length that is wanted
#' (see [find_cusum_h()] and [cusum_arl()]).
#'
#' For each signal the estimated **start of the shift** is the observation after the last time the cumulative sum was 0.
#'
#' If `target` or `sigma` is not given, the mean and the standard deviation of `x` are used. This is only appropriate when
#' `x` is a Phase I baseline; for monitoring (Phase II) supply the established values. For subgroup means with subgroup size
#' \eqn{n}, `sigma` is the standard deviation of the **mean**, \eqn{\sigma/\sqrt n}{sigma/sqrt(n)}.
#'
#' @param x Numeric vector of results (or subgroup means) in time order, with no missing values.
#' @param target In-control mean. Default: the mean of `x`.
#' @param sigma In-control standard deviation of the plotted values. Default: the standard deviation of `x`.
#' @param k Reference value in standard deviations. Default 0.5 (designed for a shift of 1 SD).
#' @param h Decision interval in standard deviations. Default 5.
#' @param head_start Starting value of both cumulative sums (a "fast initial response"), between 0 and `h`. Default 0.
#' @param reset If `TRUE`, both sums return to `head_start` after every signal, so that later signals are also found. Default `FALSE`.
#' @return An object of class `cusum_chart`: a list with `statistics` (a tibble with the columns `index`, `x`, `z`, `cplus`,
#'   `cminus`, and `signal`, which is `"up"`, `"down"`, or `""`), `signals` (a tibble with the `index`, `direction`, and the
#'   estimated `start` of each signal), `first_signal` (the index of the first signal, or `NA`), `k`, `h`, `target`,
#'   `sigma`, and `head_start`.
#' @seealso [ewma_chart()], [cusum_arl()], [find_cusum_h()], [cusum_poisson()]
#' @examples
#' set.seed(1)
#' x <- c(rnorm(30), rnorm(30, mean = 1))     # a shift of 1 SD after observation 30
#' cusum_chart(x, target = 0, sigma = 1)
#' autoplot(cusum_chart(x, target = 0, sigma = 1))
#' @export
cusum_chart <- function(x, target = NULL, sigma = NULL, k = 0.5, h = 5, head_start = 0, reset = FALSE) {
  if (!is.numeric(x) || anyNA(x) || length(x) < 2) stop("`x` must be a numeric vector with at least 2 values and no missing values.", call. = FALSE)
  if (is.null(target)) target <- mean(x) else check_finite(target, "target")
  if (is.null(sigma)) sigma <- stats::sd(x) else check_number(sigma, "sigma", positive = TRUE)
  if (!is.finite(sigma) || sigma <= 0) stop("`sigma` must be positive.", call. = FALSE)
  check_number(k, "k", positive = TRUE); check_number(h, "h", positive = TRUE)
  check_finite(head_start, "head_start")
  if (head_start < 0 || head_start >= h) stop("`head_start` must be between 0 and `h`.", call. = FALSE)
  m <- length(x)
  z <- (x - target) / sigma
  cp <- cm <- numeric(m)
  sig <- character(m)
  up <- dn <- head_start
  for (t in seq_len(m)) {
    up <- max(0, up + z[t] - k); dn <- max(0, dn - z[t] - k)
    cp[t] <- up; cm[t] <- dn
    if (up > h || dn > h) {
      sig[t] <- if (up > h && dn > h) (if (up >= dn) "up" else "down") else if (up > h) "up" else "down"
      if (reset) { up <- head_start; dn <- head_start }
    }
  }
  start_of <- function(s, path) {
    zeros <- which(path[seq_len(s)] == 0)
    if (length(zeros)) max(zeros) + 1 else 1
  }
  idx <- which(sig != "")
  starts <- vapply(idx, function(s) if (sig[s] == "up") start_of(s, cp) else start_of(s, cm), numeric(1))
  signals <- tibble::tibble(index = idx, direction = sig[idx], start = as.integer(starts))
  structure(list(statistics = tibble::tibble(index = seq_len(m), x = x, z = z, cplus = cp, cminus = cm, signal = sig),
                 signals = signals, first_signal = if (length(idx)) idx[1] else NA_integer_,
                 k = k, h = h, target = target, sigma = sigma, head_start = head_start),
            class = "cusum_chart")
}

#' @export
print.cusum_chart <- function(x, ...) {
  cat("<cusum_chart> ", nrow(x$statistics), " results, target ", format(signif(x$target, 5)), ", sigma ", format(signif(x$sigma, 4)),
      ", k = ", x$k, ", h = ", x$h, if (x$head_start > 0) paste0(", head start ", x$head_start), "\n", sep = "")
  if (is.na(x$first_signal)) {
    cat("  no signal\n")
  } else {
    s <- x$signals[1, ]
    cat("  first signal at observation ", s$index, " (", s$direction, "), the shift is estimated to start at observation ", s$start,
        "\n  ", nrow(x$signals), " observation(s) above the decision interval\n", sep = "")
  }
  invisible(x)
}

#' @export
as.data.frame.cusum_chart <- function(x, ...) as.data.frame(x$statistics)

#' Draw a CUSUM chart with ggplot2
#'
#' Plots the upper cumulative sum above zero and the lower cumulative sum below zero (as a negative value), with the decision
#' interval \eqn{\pm h} (dashed) and the signals in red.
#'
#' @param object A [cusum_chart()] object.
#' @param title,subtitle Plot title and subtitle.
#' @param xlab,ylab Axis labels.
#' @param ... Not used.
#' @return A ggplot object.
#' @export
autoplot.cusum_chart <- function(object, title = "CUSUM chart", subtitle = NULL, xlab = "Observation (time order)",
                                 ylab = "Cumulative sum (standard deviations)", ...) {
  d <- object$statistics
  d$value <- ifelse(d$cplus >= d$cminus, d$cplus, -d$cminus)
  d$status <- ifelse(d$signal == "", "No signal", "Beyond h")
  m <- nrow(d); h <- object$h
  by <- max(1, 5 * ceiling(m / 50))
  breaks <- unique(c(1, if (m >= by) seq(by, m, by = by)))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$index)) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_hline(yintercept = c(-h, h), linetype = "dashed") +
    ggplot2::geom_line(ggplot2::aes(y = .data$cplus), color = "grey60") +
    ggplot2::geom_line(ggplot2::aes(y = -.data$cminus), color = "grey60") +
    ggplot2::geom_point(ggplot2::aes(y = .data$value, color = .data$status), size = 2) +
    ggplot2::annotate("text", x = m + 0.7, y = c(h, -h), label = sprintf(c("h = %g", "-h = %g"), c(h, -h)), hjust = 0, size = 3.6) +
    ggplot2::scale_color_manual(values = c(`No signal` = "black", `Beyond h` = "red"), breaks = "Beyond h", name = NULL) +
    ggplot2::scale_x_continuous(breaks = breaks) +
    ggplot2::coord_cartesian(xlim = c(0.5, m + 0.5), clip = "off", expand = FALSE) +
    ggplot2::labs(x = xlab, y = ylab, title = title, subtitle = subtitle) +
    ggplot2::theme_bw() + ggplot2::theme(plot.margin = ggplot2::margin(6, 55, 6, 6), legend.position = "bottom")
}

#' Exponentially weighted moving average (EWMA) chart
#'
#' Computes the **EWMA** of a series of results in time order,
#' \deqn{Z_t = \lambda x_t + (1 - \lambda)Z_{t-1}, \qquad Z_0 = \mu_0,}{Z_t = lambda x_t + (1 - lambda) Z_(t-1), Z_0 = mu_0,}
#' and the control limits
#' \deqn{\mu_0 \pm L\sigma\sqrt{\frac{\lambda}{2 - \lambda}\left[1 - (1 - \lambda)^{2t}\right]}.}{mu_0 +/- L sigma sqrt(lambda/(2 - lambda) [1 - (1 - lambda)^(2t)]).}
#' The limits are narrower at the start and approach the **asymptotic** limits \eqn{\mu_0 \pm L\sigma\sqrt{\lambda/(2-\lambda)}}{mu_0 +/- L sigma sqrt(lambda/(2 - lambda))}; set
#' `limits = "asymptotic"` to use those from the first point (the ARL functions [ewma_arl()] and [sim_arl()] use them). The weight
#' \eqn{\lambda}{lambda} (between 0 and 1) controls the memory: a small \eqn{\lambda}{lambda} averages over many past results and detects small
#' shifts best, and \eqn{\lambda = 1}{lambda = 1} is a Shewhart chart.
#'
#' If `target` or `sigma` is not given, the mean and the standard deviation of `x` are used (appropriate only for a baseline).
#'
#' @param x Numeric vector of results (or subgroup means) in time order, with no missing values.
#' @param target In-control mean. Default: the mean of `x`.
#' @param sigma In-control standard deviation of the plotted values. Default: the standard deviation of `x`.
#' @param lambda Weight of the latest result, between 0 and 1. Default 0.2.
#' @param L Width of the limits in standard deviations of the EWMA. Default 3.
#' @param limits `"exact"` (default, time-varying) or `"asymptotic"`.
#' @return An object of class `ewma_chart`: a list with `statistics` (a tibble with `index`, `x`, `ewma`, `lcl`, `ucl`, and `signal`),
#'   `first_signal`, `lambda`, `L`, `target`, `sigma`, and `limits`.
#' @seealso [cusum_chart()], [ewma_arl()], [find_ewma_L()]
#' @examples
#' set.seed(1)
#' x <- c(rnorm(30), rnorm(30, mean = 1))
#' ewma_chart(x, target = 0, sigma = 1)
#' autoplot(ewma_chart(x, target = 0, sigma = 1))
#' @export
ewma_chart <- function(x, target = NULL, sigma = NULL, lambda = 0.2, L = 3, limits = c("exact", "asymptotic")) {
  limits <- match.arg(limits)
  if (!is.numeric(x) || anyNA(x) || length(x) < 2) stop("`x` must be a numeric vector with at least 2 values and no missing values.", call. = FALSE)
  if (is.null(target)) target <- mean(x) else check_finite(target, "target")
  if (is.null(sigma)) sigma <- stats::sd(x) else check_number(sigma, "sigma", positive = TRUE)
  if (!is.finite(sigma) || sigma <= 0) stop("`sigma` must be positive.", call. = FALSE)
  if (!(is.numeric(lambda) && length(lambda) == 1 && is.finite(lambda) && lambda > 0 && lambda <= 1)) {
    stop("`lambda` must be a single number in (0, 1].", call. = FALSE)
  }
  check_number(L, "L", positive = TRUE)
  m <- length(x)
  z <- numeric(m); prev <- target
  for (t in seq_len(m)) { prev <- lambda * x[t] + (1 - lambda) * prev; z[t] <- prev }
  t <- seq_len(m)
  width <- L * sigma * sqrt(lambda / (2 - lambda) * (if (limits == "exact") 1 - (1 - lambda)^(2 * t) else 1))
  sig <- z > target + width | z < target - width
  structure(list(statistics = tibble::tibble(index = t, x = x, ewma = z, lcl = target - width, ucl = target + width, signal = sig),
                 first_signal = if (any(sig)) which(sig)[1] else NA_integer_, lambda = lambda, L = L, target = target, sigma = sigma,
                 limits = limits),
            class = "ewma_chart")
}

#' @export
print.ewma_chart <- function(x, ...) {
  cat("<ewma_chart> ", nrow(x$statistics), " results, target ", format(signif(x$target, 5)), ", sigma ", format(signif(x$sigma, 4)),
      ", lambda = ", x$lambda, ", L = ", x$L, " (", x$limits, " limits)\n", sep = "")
  if (is.na(x$first_signal)) cat("  no signal\n") else
    cat("  first signal at observation ", x$first_signal, " (", sum(x$statistics$signal), " point(s) beyond the limits)\n", sep = "")
  invisible(x)
}

#' @export
as.data.frame.ewma_chart <- function(x, ...) as.data.frame(x$statistics)

#' Draw an EWMA chart with ggplot2
#'
#' Plots the EWMA, the target, the limits (dashed), and the points beyond the limits in red. The individual results are shown as
#' light points.
#'
#' @param object An [ewma_chart()] object.
#' @param title,subtitle Plot title and subtitle.
#' @param xlab,ylab Axis labels.
#' @param show_data Draw the individual results? Default `TRUE`.
#' @param ... Not used.
#' @return A ggplot object.
#' @export
autoplot.ewma_chart <- function(object, title = "EWMA chart", subtitle = NULL, xlab = "Observation (time order)",
                                ylab = "EWMA", show_data = TRUE, ...) {
  d <- object$statistics
  d$status <- ifelse(d$signal, "Beyond limits", "No signal")
  m <- nrow(d)
  by <- max(1, 5 * ceiling(m / 50))
  breaks <- unique(c(1, if (m >= by) seq(by, m, by = by)))
  g <- ggplot2::ggplot(d, ggplot2::aes(x = .data$index)) +
    ggplot2::geom_hline(yintercept = object$target)
  if (show_data) g <- g + ggplot2::geom_point(ggplot2::aes(y = .data$x), color = "grey75", size = 1.2)
  g +
    ggplot2::geom_step(ggplot2::aes(y = .data$lcl), linetype = "dashed", direction = "mid") +
    ggplot2::geom_step(ggplot2::aes(y = .data$ucl), linetype = "dashed", direction = "mid") +
    ggplot2::geom_line(ggplot2::aes(y = .data$ewma), color = "grey40") +
    ggplot2::geom_point(ggplot2::aes(y = .data$ewma, color = .data$status), size = 2) +
    ggplot2::scale_color_manual(values = c(`No signal` = "black", `Beyond limits` = "red"), breaks = "Beyond limits", name = NULL) +
    ggplot2::scale_x_continuous(breaks = breaks) +
    ggplot2::labs(x = xlab, y = ylab, title = title, subtitle = subtitle) +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
}
