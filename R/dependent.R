#' Autocorrelation function plot
#'
#' Plots the sample **autocorrelation function** (ACF) of a series in time order: the correlation between observations \eqn{k}
#' steps apart, for \eqn{k = 1, 2, \dots}, with approximate limits \eqn{\pm z/\sqrt{n}} that the autocorrelations of an
#' **independent** series stay inside (about 5% of them are outside by chance). Autocorrelations that stay large over many lags, or that
#' decay slowly, mean that the observations are dependent. For an AR(1) series with parameter \eqn{\phi} the autocorrelation at lag
#' \eqn{k} is \eqn{\phi^k}.
#'
#' @param x A numeric series in time order, with no missing values.
#' @param lag_max Largest lag. Default 30.
#' @param conf Confidence level of the limits. Default 0.95.
#' @param title,subtitle Plot title and subtitle.
#' @return A ggplot object; the autocorrelations are in its `data` (columns `lag` and `acf`).
#' @seealso [ar1_chart()], [sim_ar1()]
#' @examples
#' acf_plot(fridge_temp$temp_c, lag_max = 40)
#' @export
acf_plot <- function(x, lag_max = 30, conf = 0.95, title = "Autocorrelation function", subtitle = NULL) {
  if (!is.numeric(x) || anyNA(x) || length(x) < 5) stop("`x` must be a numeric series of at least 5 values with no missing values.", call. = FALSE)
  check_count(lag_max, "lag_max", minimum = 1)
  check_probability(conf, "conf")
  lag_max <- min(lag_max, length(x) - 2)
  a <- stats::acf(x, lag.max = lag_max, plot = FALSE)
  d <- tibble::tibble(lag = as.numeric(a$lag)[-1], acf = as.numeric(a$acf)[-1])
  lim <- stats::qnorm((1 + conf) / 2) / sqrt(length(x))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$lag, y = .data$acf)) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_hline(yintercept = c(-lim, lim), linetype = "dashed", color = "#4c78a8") +
    ggplot2::geom_segment(ggplot2::aes(xend = .data$lag, yend = 0), linewidth = 1.1) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::scale_y_continuous(limits = c(min(-0.2, min(d$acf) - 0.05), 1)) +
    ggplot2::labs(x = "Lag (number of observations apart)", y = "Autocorrelation", title = title, subtitle = subtitle) +
    ggplot2::theme_bw() + ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

#' Simulate an AR(1) series
#'
#' Simulates the first-order autoregressive series \eqn{x_t - \mu = \phi(x_{t-1} - \mu) + e_t} with independent normal innovations
#' \eqn{e_t}, started in its stationary distribution. The **marginal standard deviation** of the series is `sigma` (so the innovations have
#' standard deviation \eqn{\sigma\sqrt{1 - \phi^2}}), and the autocorrelation at lag \eqn{k} is \eqn{\phi^k}. A step shift of the mean can be added.
#'
#' @param n Length of the series.
#' @param phi Autoregressive parameter, between -1 and 1 (the lag-1 autocorrelation).
#' @param sigma Marginal standard deviation. Default 1.
#' @param mean Mean of the series. Default 0.
#' @param shift,change A step shift of size `shift` starting at observation `change` (both must be given). Default no shift.
#' @return A numeric vector of length `n`. Use [set.seed()] before the call for a reproducible result.
#' @seealso [sim_dependent()], [ar1_chart()]
#' @examples
#' set.seed(1)
#' x <- sim_ar1(200, phi = 0.9, sigma = 0.5, mean = 5)
#' acf_plot(x)
#' @export
sim_ar1 <- function(n, phi, sigma = 1, mean = 0, shift = 0, change = NULL) {
  check_count(n, "n", minimum = 2)
  if (!(is.numeric(phi) && length(phi) == 1 && is.finite(phi) && abs(phi) < 1)) stop("`phi` must be a single number between -1 and 1.", call. = FALSE)
  check_number(sigma, "sigma", positive = TRUE); check_finite(mean, "mean"); check_finite(shift, "shift")
  if (shift != 0) {
    if (is.null(change)) stop("Give `change`, the observation at which the shift starts.", call. = FALSE)
    check_count(change, "change", minimum = 1)
  }
  e <- stats::rnorm(n, 0, sigma * sqrt(1 - phi^2))
  x <- numeric(n); x[1] <- stats::rnorm(1, 0, sigma)
  for (t in seq_len(n)[-1]) x[t] <- phi * x[t - 1] + e[t]
  x <- x + mean
  if (shift != 0) x[seq(change, n)] <- x[seq(change, n)] + shift
  x
}

#' The false alarm probability of an individuals chart for AR(1) data
#'
#' The moving range estimates the standard deviation of the **difference** of consecutive values, \eqn{\sigma\sqrt{2(1 - \rho)}}. So for a series
#' with lag-1 autocorrelation \eqn{\rho} the estimate \eqn{\overline{MR}/1.128} is about \eqn{\sigma\sqrt{1 - \rho}}, and the limits of the I chart,
#' \eqn{\bar x \pm k\overline{MR}/1.128}, are at \eqn{\pm k\sqrt{1 - \rho}} marginal standard deviations. This function gives the probability that a
#' value of a normal stationary series is beyond these limits (with the mean known): \eqn{2\{1 - \Phi(k\sqrt{1 - \rho})\}}.
#'
#' @param rho Lag-1 autocorrelation (the AR(1) parameter), between -1 and 1.
#' @param k Number of standard deviations of the limits. Default 3.
#' @return The probability that a point is beyond the limits.
#' @examples
#' i_chart_alarm_rate(c(0, 0.5, 0.9))
#' @export
i_chart_alarm_rate <- function(rho, k = 3) {
  if (!(is.numeric(rho) && all(is.finite(rho)) && all(abs(rho) < 1))) stop("`rho` must be numbers between -1 and 1.", call. = FALSE)
  check_number(k, "k", positive = TRUE)
  2 * stats::pnorm(-k * sqrt(1 - rho))
}

#' Residual chart for an AR(1) series
#'
#' For a series whose consecutive values are **dependent**, charting the values gives limits that are too narrow and many false alarms. If the series
#' follows the AR(1) model \eqn{x_t - \mu = \phi(x_{t-1} - \mu) + e_t}, the **one-step-ahead residuals** \eqn{e_t = (x_t - \mu) - \phi(x_{t-1} - \mu)} are
#' independent with a constant standard deviation, so an individuals chart of the residuals (a **special-cause chart**) has the usual false alarm rate. A
#' residual beyond \eqn{\pm k\hat\sigma_e} says that the value is **unexpected given the previous one**. The parameters \eqn{\mu} and \eqn{\phi} are
#' estimated with [stats::arima()] unless they are given (for Phase II monitoring); \eqn{\sigma_e} is estimated from the average moving range of the
#' residuals (`sigma_from = "mr"`, as for the I chart) or from their root mean square (`"sd"`) unless it is given.
#'
#' A **sustained** shift of the mean by \eqn{\delta} changes the residuals by only \eqn{(1 - \phi)\delta} after the first point, so a residual chart is
#' much more sensitive to sudden departures than to slow drifts when \eqn{\phi} is large; a CUSUM or EWMA of the residuals, or a fixed specification limit,
#' covers slow drifts.
#'
#' @param x A numeric series in time order, with no missing values.
#' @param phi,mean,sigma_e Known AR(1) parameter, mean, and innovation standard deviation. Default: estimated.
#' @param k Number of standard deviations of the limits. Default 3.
#' @param sigma_from `"mr"` (default) or `"sd"`: how \eqn{\sigma_e} is estimated from the residuals.
#' @return An object of class `ar1_chart`: a list with `statistics` (a tibble with `index`, `x`, `fitted` (the one-step-ahead prediction), `residual`,
#'   `lcl`, `ucl` (limits for the residual), and `signal`), `signals` (the indices of the points beyond the limits), `phi`, `mean`, `sigma_e`, `sigma_x`
#'   (the marginal standard deviation \eqn{\sigma_e/\sqrt{1 - \phi^2}}), and `k`.
#' @seealso [acf_plot()], [qc_chart()] for the I chart, [sim_dependent()]
#' @examples
#' ar1_chart(fridge_temp$temp_c, k = 4)
#' autoplot(ar1_chart(fridge_temp$temp_c, k = 4))
#' @export
ar1_chart <- function(x, phi = NULL, mean = NULL, sigma_e = NULL, k = 3, sigma_from = c("mr", "sd")) {
  sigma_from <- match.arg(sigma_from)
  min_n <- if (is.null(phi) || is.null(mean)) 10 else 3
  if (!is.numeric(x) || anyNA(x) || length(x) < min_n) stop("`x` must be a numeric series of at least ", min_n, " values with no missing values.", call. = FALSE)
  check_number(k, "k", positive = TRUE)
  if (is.null(phi) || is.null(mean)) {
    fit <- tryCatch(stats::arima(x, order = c(1, 0, 0), include.mean = TRUE), error = function(e) NULL)
    if (is.null(fit)) stop("The AR(1) model could not be fitted; check that the series is stationary.", call. = FALSE)
    if (is.null(phi)) phi <- unname(stats::coef(fit)[["ar1"]])
    if (is.null(mean)) mean <- unname(stats::coef(fit)[["intercept"]])
  }
  if (!(is.numeric(phi) && length(phi) == 1 && is.finite(phi) && abs(phi) < 1)) stop("`phi` must be a single number between -1 and 1.", call. = FALSE)
  check_finite(mean, "mean")
  n <- length(x)
  fitted <- c(NA_real_, mean + phi * (x[-n] - mean))
  res <- x - fitted
  if (is.null(sigma_e)) {
    r <- res[-1]
    sigma_e <- if (sigma_from == "mr") mean(abs(diff(r))) / 1.128 else sqrt(mean(r^2))
  } else check_number(sigma_e, "sigma_e", positive = TRUE)
  sig <- !is.na(res) & abs(res) > k * sigma_e
  stat <- tibble::tibble(index = seq_len(n), x = x, fitted = fitted, residual = res, lcl = -k * sigma_e, ucl = k * sigma_e, signal = sig)
  structure(list(statistics = stat, signals = which(sig), phi = phi, mean = mean, sigma_e = sigma_e, sigma_x = sigma_e / sqrt(1 - phi^2), k = k),
            class = "ar1_chart")
}

#' @export
print.ar1_chart <- function(x, ...) {
  cat("<ar1_chart> AR(1) fit: phi = ", format(round(x$phi, 3)), ", mean = ", format(signif(x$mean, 5)), ", sigma_e = ", format(signif(x$sigma_e, 3)),
      " (marginal SD ", format(signif(x$sigma_x, 3)), ")\n", sep = "")
  if (length(x$signals) == 0) cat("  no residual beyond the ", x$k, "-sigma limits\n", sep = "") else
    cat("  ", length(x$signals), " residual(s) beyond the ", x$k, "-sigma limits: observation(s) ", paste(utils::head(x$signals, 10), collapse = ", "),
        if (length(x$signals) > 10) ", ...", "\n", sep = "")
  invisible(x)
}

#' @export
as.data.frame.ar1_chart <- function(x, ...) as.data.frame(x$statistics)

#' Draw an AR(1) residual chart with ggplot2
#'
#' `which = "residual"` draws the individuals chart of the residuals. `which = "data"` draws the series with the **prediction band**
#' (the one-step-ahead prediction plus or minus the limits): a value outside the band is a residual beyond the limits.
#'
#' @param object An [ar1_chart()] object.
#' @param which `"residual"` (default) or `"data"`.
#' @param title,subtitle Plot title and subtitle.
#' @param xlab,ylab Axis labels.
#' @param ... Not used.
#' @return A ggplot object.
#' @export
autoplot.ar1_chart <- function(object, which = c("residual", "data"), title = NULL, subtitle = NULL, xlab = "Observation (time order)", ylab = NULL, ...) {
  which <- match.arg(which)
  d <- object$statistics
  d$status <- ifelse(d$signal, "Beyond limits", "No signal")
  m <- nrow(d); by <- max(1, 5 * ceiling(m / 50))
  g <- ggplot2::ggplot(d, ggplot2::aes(x = .data$index))
  if (which == "residual") {
    if (is.null(title)) title <- "Residual chart of the AR(1) fit"
    if (is.null(ylab)) ylab <- "Residual"
    g <- g + ggplot2::geom_hline(yintercept = 0) + ggplot2::geom_hline(yintercept = c(object$statistics$lcl[1], object$statistics$ucl[1]), linetype = "dashed") +
      ggplot2::geom_line(ggplot2::aes(y = .data$residual), color = "grey60", na.rm = TRUE) +
      ggplot2::geom_point(ggplot2::aes(y = .data$residual, color = .data$status), size = 1.6, na.rm = TRUE)
  } else {
    if (is.null(title)) title <- "Series with the one-step-ahead prediction band"
    if (is.null(ylab)) ylab <- "Value"
    g <- g + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$fitted + .data$lcl, ymax = .data$fitted + .data$ucl), fill = "#4c78a8", alpha = 0.25, na.rm = TRUE) +
      ggplot2::geom_line(ggplot2::aes(y = .data$x), color = "grey40") +
      ggplot2::geom_point(ggplot2::aes(y = .data$x, color = .data$status), size = 1.4)
  }
  g + ggplot2::scale_color_manual(values = c(`No signal` = "black", `Beyond limits` = "red"), breaks = "Beyond limits", name = NULL) +
    ggplot2::labs(x = xlab, y = ylab, title = title, subtitle = subtitle) +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
}

#' Simulate the false alarm rate of control charts for AR(1) data
#'
#' Simulates `reps` in-control AR(1) series of length `n` (mean 0, marginal standard deviation 1) and builds, from each series, three charts, counting
#' the points beyond their \eqn{\pm k}-sigma limits:
#' * `"I"`: the individuals chart with the standard estimate \eqn{\hat\sigma = \overline{MR}/1.128}, which is biased low for positive \eqn{\phi};
#' * `"I_sd"`: the individuals chart with \eqn{\hat\sigma} the sample standard deviation of the series (the limits are right for the **marginal**
#'   distribution, but the alarms come in clusters);
#' * `"residual"`: the residual chart of an AR(1) model fitted to the series, as in [ar1_chart()] but with the Yule-Walker (lag-1
#'   sample autocorrelation) estimate of \eqn{\phi} instead of [stats::arima()], and \eqn{\sigma_e} from the average moving range of the
#'   residuals.
#'
#' @param phi AR(1) parameter (the lag-1 autocorrelation), a single number between -1 and 1.
#' @param n Length of each series. Default 200.
#' @param reps Number of series. Default 1000.
#' @param k Number of standard deviations. Default 3.
#' @param charts Which charts to evaluate. Default all three.
#' @return A tibble with one row per series and chart, with `rep` (the series), `chart`, `out` (the number of points beyond the limits),
#'   `share` (`out / n`, or `out / (n - 1)` for the residual chart, which has no residual for the first point), and `any` (at least one
#'   point beyond the limits).
#' @seealso [i_chart_alarm_rate()], [sim_ar1()]
#' @examples
#' set.seed(1)
#' s <- sim_dependent(phi = 0.9, reps = 200)
#' tapply(s$share, s$chart, mean)
#' @export
sim_dependent <- function(phi, n = 200, reps = 1000, k = 3, charts = c("I", "I_sd", "residual")) {
  charts <- match.arg(charts, several.ok = TRUE)
  check_count(n, "n", minimum = 10); check_count(reps, "reps", minimum = 1); check_number(k, "k", positive = TRUE)
  out <- vector("list", reps * length(charts)); j <- 0
  for (r in seq_len(reps)) {
    x <- sim_ar1(n, phi)
    xc <- x - mean(x)
    for (ch in charts) {
      o <- switch(ch,
        I = sum(abs(xc) > k * mean(abs(diff(x))) / 1.128),
        I_sd = sum(abs(xc) > k * stats::sd(x)),
        residual = { ph <- sum(xc[-1] * xc[-n]) / sum(xc^2); e <- xc[-1] - ph * xc[-n]; sum(abs(e) > k * mean(abs(diff(e))) / 1.128) })
      j <- j + 1
      out[[j]] <- tibble::tibble(rep = r, chart = ch, out = o, share = o / (if (ch == "residual") n - 1 else n), any = o > 0)
    }
  }
  do.call(rbind, out)
}
