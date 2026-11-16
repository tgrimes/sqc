#' Control charts for counts and proportions: p, np, c, and u charts
#'
#' Computes a Shewhart chart for **attribute data**. For \eqn{m} samples with counts \eqn{x_i} (and, for the p, np, and u
#' charts, sample sizes or exposures \eqn{n_i}), the chart type is chosen by what is counted:
#'
#' | `type` | counted | plotted | center | standard deviation of the count |
#' |:--|:--|:--|:--|:--|
#' | `"p"` | nonconforming units out of \eqn{n_i} (varying) | \eqn{x_i/n_i} | \eqn{\bar p = \sum x_i/\sum n_i} | \eqn{\sqrt{n_i\bar p(1-\bar p)}} |
#' | `"np"` | nonconforming units out of a constant \eqn{n} | \eqn{x_i} | \eqn{n\bar p} | \eqn{\sqrt{n\bar p(1-\bar p)}} |
#' | `"c"` | events in a constant area of opportunity | \eqn{x_i} | \eqn{\bar c = \bar x} | \eqn{\sqrt{\bar c}} |
#' | `"u"` | events in a varying exposure \eqn{n_i} (line-days, patient-days) | \eqn{x_i/n_i} | \eqn{\bar u = \sum x_i/\sum n_i} | \eqn{\sqrt{n_i\bar u}} |
#'
#' Because the variance of a count depends on its mean and on the sample size, the limits of the p and u charts **vary with
#' \eqn{n_i}**. With `limits = "normal"` the limits are the center line plus or minus `k` standard deviations (truncated at 0
#' and, for the p chart, at 1). With `limits = "exact"` they are **probability limits** from the binomial or Poisson
#' distribution: a point signals if the count is above the \eqn{1 - \Phi(-k)} quantile or below the \eqn{\Phi(-k)} quantile, so
#' that each side has a false alarm probability of at most \eqn{\Phi(-k)} (0.00135 for `k = 3`), which the normal limits do not
#' achieve for small expected counts.
#'
#' With `adjust = "laney"` the standard deviation is multiplied by \eqn{\sigma_z}, the standard deviation of the
#' standardized values \eqn{z_i} estimated from their average moving range (\eqn{\overline{MR}/1.128}), which widens the limits
#' when the data are **overdispersed** (the Laney p' and u' charts).
#'
#' @param count Non-negative whole numbers: the counts, in time order.
#' @param n Sample sizes (p, np) or exposures (u), one value or one per count. Not used for the c chart.
#' @param type `"p"`, `"np"`, `"c"`, or `"u"`.
#' @param center Optional known center: the proportion (p, np), the mean count (c), or the rate (u). Default: estimated from
#'   the data (appropriate for a Phase I baseline).
#' @param k Number of standard deviations (or the corresponding tail probability, for exact limits). Default 3.
#' @param limits `"normal"` (default) or `"exact"`.
#' @param adjust `"none"` (default) or `"laney"`.
#' @return An object of class `count_chart`: a list with `statistics` (a tibble with `index`, `count`, `n`, `value`, `center`,
#'   `lcl`, `ucl`, and `signal`), `signals` (the points beyond the limits, with the side), `type`, `center`, `k`, `limits`,
#'   `adjust`, and `sigma_z`.
#' @seealso [funnel_plot()] for comparing many institutions, [cusum_poisson()] for counts in time order.
#' @examples
#' count_chart(hand_hygiene$compliant, hand_hygiene$audits, type = "p")
#' u <- count_chart(clabsi_unit$infections, clabsi_unit$line_days / 1000, type = "u",
#'                  limits = "exact")
#' autoplot(u, ylab = "Infections per 1000 line-days")
#' @export
count_chart <- function(count, n = NULL, type = c("p", "np", "c", "u"), center = NULL, k = 3,
                        limits = c("normal", "exact"), adjust = c("none", "laney")) {
  type <- match.arg(type); limits <- match.arg(limits); adjust <- match.arg(adjust)
  if (!is.numeric(count) || anyNA(count) || any(count < 0) || any(count != round(count)) || length(count) < 2) {
    stop("`count` must be at least 2 non-negative whole numbers with no missing values.", call. = FALSE)
  }
  check_number(k, "k", positive = TRUE)
  m <- length(count)
  if (type == "c") {
    n <- rep(1, m)
  } else {
    if (is.null(n)) stop("`n` (the sample sizes or exposures) is required for a ", type, " chart.", call. = FALSE)
    if (!is.numeric(n) || anyNA(n) || any(n <= 0) || !(length(n) %in% c(1, m))) {
      stop("`n` must be positive numbers, one value or one per count.", call. = FALSE)
    }
    if (length(n) == 1) n <- rep(n, m)
  }
  if (type %in% c("p", "np") && any(count > n)) stop("`count` cannot exceed `n` for a ", type, " chart.", call. = FALSE)
  if (type == "np" && length(unique(n)) > 1) stop("An np chart needs a constant `n`; use a p chart for varying sample sizes.", call. = FALSE)
  if (adjust == "laney" && limits == "exact") stop("The Laney adjustment is defined for normal limits.", call. = FALSE)

  binom <- type %in% c("p", "np")
  if (is.null(center)) {
    center <- if (type == "c") mean(count) else sum(count) / sum(n)
  } else {
    check_number(center, "center", positive = TRUE)
    if (binom && center >= 1) stop("`center` must be a proportion below 1 for this chart.", call. = FALSE)
  }
  if (type == "c" && center <= 0) stop("The center must be positive.", call. = FALSE)
  if (center <= 0) stop("The estimated center is 0; the limits are not defined.", call. = FALSE)

  mu <- switch(type, p = n * center, np = n * center, c = rep(center, m), u = n * center)  # expected counts
  sd_count <- if (binom) sqrt(n * center * (1 - center)) else sqrt(mu)
  per_n <- type %in% c("p", "u")                                                         # plotted value is count / n
  value <- if (per_n) count / n else count
  z <- (count - mu) / sd_count
  sigma_z <- if (adjust == "laney" && m >= 3) mean(abs(diff(z))) / 1.128 else 1
  if (limits == "normal") {
    lo <- pmax(0, mu - k * sigma_z * sd_count); hi <- mu + k * sigma_z * sd_count
    if (binom) hi <- pmin(hi, n)
  } else {
    a <- stats::pnorm(-k)
    if (binom) { lo <- stats::qbinom(a, n, center); hi <- stats::qbinom(1 - a, n, center) }
    else { lo <- stats::qpois(a, mu); hi <- stats::qpois(1 - a, mu) }
  }
  side <- ifelse(count > hi, "high", ifelse(count < lo, "low", ""))
  cl <- if (per_n) rep(center, m) else mu
  stat <- tibble::tibble(index = seq_len(m), count = count, n = if (type == "c") NA_real_ else n, value = value, center = cl,
                         lcl = if (per_n) lo / n else lo, ucl = if (per_n) hi / n else hi, signal = side != "")
  structure(list(statistics = stat, signals = tibble::tibble(index = which(side != ""), value = value[side != ""], side = side[side != ""]),
                 type = type, center = center, k = k, limits = limits, adjust = adjust, sigma_z = sigma_z),
            class = "count_chart")
}

#' @export
print.count_chart <- function(x, ...) {
  nm <- c(p = "p", np = "np", c = "c", u = "u")[[x$type]]
  cat("<count_chart> ", nm, " chart, ", nrow(x$statistics), " samples, center ", format(signif(x$center, 4)), ", ",
      if (x$limits == "exact") paste0("exact ", x$k, "-sigma probability limits") else paste0(x$k, "-sigma limits"),
      if (x$adjust == "laney") paste0(" (Laney, sigma_z = ", format(round(x$sigma_z, 2)), ")"), "\n", sep = "")
  if (nrow(x$signals) == 0) cat("  no signals\n") else
    cat("  ", nrow(x$signals), " signal(s): sample(s) ", paste(x$signals$index, collapse = ", "), "\n", sep = "")
  invisible(x)
}

#' @export
as.data.frame.count_chart <- function(x, ...) as.data.frame(x$statistics)

#' Draw a count chart with ggplot2
#'
#' Plots the rate or count, the center line, and the limits (as steps, since they vary with the sample size), with the signals in red.
#'
#' @param object A [count_chart()] object.
#' @param title,subtitle Plot title and subtitle.
#' @param xlab,ylab Axis labels.
#' @param ... Not used.
#' @return A ggplot object.
#' @export
autoplot.count_chart <- function(object, title = NULL, subtitle = NULL, xlab = "Sample", ylab = NULL, ...) {
  d <- object$statistics
  d$status <- ifelse(d$signal, "Beyond limits", "No signal")
  if (is.null(title)) title <- paste0(c(p = "p", np = "np", c = "c", u = "u")[[object$type]], " chart")
  if (is.null(ylab)) ylab <- switch(object$type, p = "Proportion", np = "Count", c = "Count", u = "Rate")
  m <- nrow(d)
  by <- max(1, 5 * ceiling(m / 50))
  breaks <- unique(c(1, if (m >= by) seq(by, m, by = by)))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$index)) +
    ggplot2::geom_step(ggplot2::aes(y = .data$ucl), direction = "mid", linetype = "dashed") +
    ggplot2::geom_step(ggplot2::aes(y = .data$lcl), direction = "mid", linetype = "dashed") +
    ggplot2::geom_line(ggplot2::aes(y = .data$center)) +
    ggplot2::geom_line(ggplot2::aes(y = .data$value), color = "grey60") +
    ggplot2::geom_point(ggplot2::aes(y = .data$value, color = .data$status), size = 2) +
    ggplot2::scale_color_manual(values = c(`No signal` = "black", `Beyond limits` = "red"), breaks = "Beyond limits", name = NULL) +
    ggplot2::scale_x_continuous(breaks = breaks) +
    ggplot2::labs(x = xlab, y = ylab, title = title, subtitle = subtitle) +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
}
