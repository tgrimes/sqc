#' Likelihood-ratio CUSUMs for counts and for risk-adjusted outcomes
#'
#' A CUSUM can be built from the **log-likelihood ratio** of two hypotheses about each new observation: that the process
#' is as expected, and that it has changed by a stated amount. Adding the log-likelihood ratios \eqn{W_t} and never letting the
#' sum go below 0,
#' \deqn{S_t = \max(0,\; S_{t-1} + W_t),}{S_t = max(0, S_(t-1) + W_t),}
#' gives a chart that signals when \eqn{S_t > h}. The two charts here differ only in the weights \eqn{W_t}:
#' * `cusum_poisson()` is the **observed-minus-expected** CUSUM for counts (for example deaths against the deaths expected for the
#'   age and sex of the patients). With observed counts \eqn{O_t} and expected counts \eqn{E_t}, and a relative risk `rr` to
#'   detect, \eqn{W_t = O_t\log(rr) - E_t(rr - 1)}{W_t = O_t log(rr) - E_t (rr - 1)}. With `rr > 1` the chart detects an increase, with
#'   `rr < 1` a decrease.
#' * `cusum_risk()` is the **risk-adjusted Bernoulli CUSUM** of Steiner et al. (2000) for binary outcomes \eqn{y_t} (1 = the event, for
#'   example death) with a predicted risk \eqn{p_t} for each patient. For an odds ratio `odds_ratio` (OR) to detect,
#'   \eqn{W_t = y_t\log(OR) - \log(1 - p_t + OR\,p_t)}{W_t = y_t log(OR) - log(1 - p_t + OR p_t)}. The predicted risks are **given**: nothing is fitted.
#'
#' The decision interval `h` can be left out to look at the statistic; choose it by simulation ([sim_cusum_poisson()],
#' [sim_cusum_risk()], [cusum_limit()]). For each signal the **estimated start of the change** is the observation after the last time
#' the sum was 0, and the **accumulated excess** is the sum of observed minus expected over the observations since then.
#'
#' @param observed Observed counts (non-negative whole numbers) in time order.
#' @param expected Expected counts under the in-control model (positive), the same length as `observed`.
#' @param outcome Binary outcomes, 0 or 1, in time order.
#' @param risk Predicted probabilities of the event (between 0 and 1), the same length as `outcome`.
#' @param rr Relative risk to detect (positive, not 1). Default 1.5.
#' @param odds_ratio Odds ratio to detect (positive, not 1). Default 2.
#' @param h Decision interval; the chart signals when the statistic exceeds it. Default `NULL` (no limit).
#' @param reset If `TRUE`, the statistic returns to 0 after each signal. Default `FALSE`.
#' @return An object of class `cusum_llr`: a list with `statistics` (a tibble with `index`, `observed`, `expected`, `weight`,
#'   `statistic`, and `signal`; for `cusum_risk()` the `observed` column holds the outcomes and `expected` the predicted risks),
#'   `signals` (a tibble with the `index`, the estimated `start`, and the accumulated `excess` of each signal), `first_signal`, `h`,
#'   `kind` (`"poisson"` or `"risk"`), and `design` (the design value as a named number, `c(rr = )` or `c(odds_ratio = )`).
#' @references Steiner, S. H., Cook, R. J., Farewell, V. T. and Treasure, T. (2000). Monitoring surgical performance using risk-adjusted
#'   cumulative sum charts. *Biostatistics*, 1, 441-452. Spiegelhalter, D., Grigg, O., Kinsman, R. and Treasure, T. (2003).
#'   Risk-adjusted sequential probability ratio tests. *International Journal for Quality in Health Care*, 15, 7-13.
#' @seealso [cusum_chart()], [sim_cusum_poisson()], [sim_cusum_risk()]
#' @examples
#' cusum_poisson(gp_deaths$observed, gp_deaths$expected, rr = 1.5, h = 4)
#' surgeon_c <- surgeon_cases[surgeon_cases$surgeon == "C", ]
#' cusum_risk(surgeon_c$death, surgeon_c$risk, odds_ratio = 2, h = 3)
#' @export
cusum_poisson <- function(observed, expected, rr = 1.5, h = NULL, reset = FALSE) {
  if (!is.numeric(observed) || anyNA(observed) || any(observed < 0) || any(observed != round(observed))) {
    stop("`observed` must be non-negative whole numbers with no missing values.", call. = FALSE)
  }
  if (!is.numeric(expected) || anyNA(expected) || any(expected <= 0) || length(expected) != length(observed)) {
    stop("`expected` must be positive numbers, one for each observed count.", call. = FALSE)
  }
  check_number(rr, "rr", positive = TRUE)
  if (rr == 1) stop("`rr` must differ from 1.", call. = FALSE)
  w <- observed * log(rr) - expected * (rr - 1)
  make_llr(w, observed, expected, h, reset, kind = "poisson", design = c(rr = rr))
}

#' @rdname cusum_poisson
#' @export
cusum_risk <- function(outcome, risk, odds_ratio = 2, h = NULL, reset = FALSE) {
  if (!is.numeric(outcome) || anyNA(outcome) || !all(outcome %in% c(0, 1))) stop("`outcome` must be 0 or 1 with no missing values.", call. = FALSE)
  if (!is.numeric(risk) || anyNA(risk) || any(risk <= 0 | risk >= 1) || length(risk) != length(outcome)) {
    stop("`risk` must be probabilities strictly between 0 and 1, one for each outcome.", call. = FALSE)
  }
  check_number(odds_ratio, "odds_ratio", positive = TRUE)
  if (odds_ratio == 1) stop("`odds_ratio` must differ from 1.", call. = FALSE)
  w <- outcome * log(odds_ratio) - log(1 - risk + odds_ratio * risk)
  make_llr(w, outcome, risk, h, reset, kind = "risk", design = c(odds_ratio = odds_ratio))
}

make_llr <- function(w, observed, expected, h, reset, kind, design) {
  if (!is.null(h)) check_number(h, "h", positive = TRUE)
  m <- length(w)
  s <- numeric(m); sig <- logical(m); prev <- 0
  for (t in seq_len(m)) {
    prev <- max(0, prev + w[t]); s[t] <- prev
    if (!is.null(h) && prev > h) { sig[t] <- TRUE; if (reset) prev <- 0 }
  }
  idx <- which(sig)
  start <- vapply(idx, function(t) { z <- which(s[seq_len(t)] == 0); if (length(z)) max(z) + 1 else 1 }, numeric(1))
  excess <- mapply(function(a, b) sum(observed[a:b] - expected[a:b]), start, idx)
  structure(list(statistics = tibble::tibble(index = seq_len(m), observed = observed, expected = expected, weight = w, statistic = s, signal = sig),
                 signals = tibble::tibble(index = idx, start = as.integer(start), excess = as.numeric(excess)),
                 first_signal = if (length(idx)) idx[1] else NA_integer_, h = h, kind = kind, design = design),
            class = "cusum_llr")
}

#' @export
print.cusum_llr <- function(x, ...) {
  what <- if (x$kind == "poisson") "observed-minus-expected (Poisson)" else "risk-adjusted (Bernoulli)"
  cat("<cusum_llr> ", what, " CUSUM, ", nrow(x$statistics), " observations, ", names(x$design), " = ", x$design,
      if (!is.null(x$h)) paste0(", h = ", x$h), "\n", sep = "")
  cat("  maximum of the statistic: ", format(round(max(x$statistics$statistic), 2)), " at observation ", which.max(x$statistics$statistic), "\n", sep = "")
  if (is.null(x$h)) {
    cat("  no decision interval given\n")
  } else if (is.na(x$first_signal)) {
    cat("  no signal\n")
  } else {
    s <- x$signals[1, ]
    cat("  first signal at observation ", s$index, ": the change is estimated to start at observation ", s$start,
        ", with an accumulated excess of ", format(round(s$excess, 1)), " (observed minus expected)\n", sep = "")
  }
  invisible(x)
}

#' @export
as.data.frame.cusum_llr <- function(x, ...) as.data.frame(x$statistics)

#' Draw a likelihood-ratio CUSUM with ggplot2
#'
#' Plots the CUSUM statistic in time order, the decision interval `h` (dashed) if there is one, and the signals in red.
#'
#' @param object A [cusum_poisson()] or [cusum_risk()] object.
#' @param title,subtitle Plot title and subtitle.
#' @param xlab,ylab Axis labels.
#' @param ... Not used.
#' @return A ggplot object.
#' @export
autoplot.cusum_llr <- function(object, title = NULL, subtitle = NULL, xlab = "Observation (time order)", ylab = "CUSUM statistic", ...) {
  if (is.null(title)) title <- if (object$kind == "poisson") "Observed-minus-expected CUSUM" else "Risk-adjusted CUSUM"
  d <- object$statistics
  d$status <- ifelse(d$signal, "Beyond h", "No signal")
  m <- nrow(d)
  by <- max(1, 5 * ceiling(m / 50))
  breaks <- unique(c(1, if (m >= by) seq(by, m, by = by)))
  g <- ggplot2::ggplot(d, ggplot2::aes(x = .data$index, y = .data$statistic)) +
    ggplot2::geom_line(color = "grey60") +
    ggplot2::geom_point(ggplot2::aes(color = .data$status), size = 1.8)
  if (!is.null(object$h)) {
    g <- g + ggplot2::geom_hline(yintercept = object$h, linetype = "dashed") +
      ggplot2::annotate("text", x = m + 0.7, y = object$h, label = sprintf("h = %g", object$h), hjust = 0, size = 3.6)
  }
  g +
    ggplot2::scale_color_manual(values = c(`No signal` = "black", `Beyond h` = "red"), breaks = "Beyond h", name = NULL) +
    ggplot2::scale_x_continuous(breaks = breaks) +
    ggplot2::coord_cartesian(xlim = c(0.5, m + 0.5), clip = "off", expand = FALSE) +
    ggplot2::labs(x = xlab, y = ylab, title = title, subtitle = subtitle) +
    ggplot2::theme_bw() + ggplot2::theme(plot.margin = ggplot2::margin(6, 45, 6, 6), legend.position = "bottom")
}

# Run the CUSUM recursion down the columns of a matrix of weights (one row per simulated series).
llr_paths <- function(W, h) {
  reps <- nrow(W); horizon <- ncol(W)
  s <- numeric(reps); mx <- numeric(reps); rl <- rep(NA_integer_, reps)
  for (t in seq_len(horizon)) {
    s <- pmax(0, s + W[, t])
    mx <- pmax(mx, s)
    if (!is.null(h)) { hit <- is.na(rl) & s > h; rl[hit] <- t }
  }
  tibble::tibble(max_stat = mx, run_length = rl)
}

#' Simulate a likelihood-ratio CUSUM to choose its decision interval
#'
#' With the expected counts (`sim_cusum_poisson()`) or predicted risks (`sim_cusum_risk()`) of the series to be monitored, simulates
#' many series from the **true** model, runs the CUSUM of [cusum_poisson()] or [cusum_risk()] on each, and records the largest value
#' of the statistic and (if `h` is given) the observation at which it first exceeds `h`. With `true_rr = 1` (or `true_odds_ratio = 1`)
#' the series are in control, and the share of series whose maximum exceeds `h` is the **false alarm probability** over the
#' length of the series. [cusum_limit()] finds the `h` that gives a chosen false alarm probability. With a true change, the
#' run lengths show how quickly the chart signals.
#'
#' The count series are Poisson with mean `true_rr * expected`; the outcomes are Bernoulli with the predicted risks changed to the odds ratio
#' `true_odds_ratio`. Use [set.seed()] before the call for a reproducible result.
#'
#' @param expected Expected counts of the series (the length of the series is the horizon).
#' @param risk Predicted risks of the series.
#' @param rr,odds_ratio The relative risk or odds ratio that the chart is designed to detect. Defaults 1.5 and 2.
#' @param true_rr,true_odds_ratio The true relative risk or odds ratio used in the simulation. Default 1 (in control).
#' @param h Optional decision interval, to record run lengths.
#' @param reps Number of simulated series. Default 5000.
#' @return A tibble with one row per series and columns `max_stat` and `run_length` (`NA` if there was no signal or no `h`).
#' @seealso [cusum_poisson()], [cusum_risk()], [cusum_limit()]
#' @examples
#' set.seed(1)
#' s <- sim_cusum_poisson(gp_deaths$expected, rr = 1.5, reps = 1000)
#' cusum_limit(s, alpha = 0.05)       # the h with a 5% false alarm probability over the years
#' @export
sim_cusum_poisson <- function(expected, rr = 1.5, true_rr = 1, h = NULL, reps = 5000) {
  if (!is.numeric(expected) || anyNA(expected) || any(expected <= 0) || length(expected) < 1) stop("`expected` must be positive numbers.", call. = FALSE)
  check_number(rr, "rr", positive = TRUE); check_number(true_rr, "true_rr", positive = TRUE)
  if (rr == 1) stop("`rr` must differ from 1.", call. = FALSE)
  if (!is.null(h)) check_number(h, "h", positive = TRUE)
  check_count(reps, "reps", minimum = 1)
  horizon <- length(expected)
  O <- matrix(stats::rpois(reps * horizon, rep(true_rr * expected, each = reps)), reps, horizon)
  W <- O * log(rr) - rep(expected * (rr - 1), each = reps)
  llr_paths(W, h)
}

#' @rdname sim_cusum_poisson
#' @export
sim_cusum_risk <- function(risk, odds_ratio = 2, true_odds_ratio = 1, h = NULL, reps = 5000) {
  if (!is.numeric(risk) || anyNA(risk) || any(risk <= 0 | risk >= 1) || length(risk) < 1) stop("`risk` must be probabilities strictly between 0 and 1.", call. = FALSE)
  check_number(odds_ratio, "odds_ratio", positive = TRUE); check_number(true_odds_ratio, "true_odds_ratio", positive = TRUE)
  if (odds_ratio == 1) stop("`odds_ratio` must differ from 1.", call. = FALSE)
  if (!is.null(h)) check_number(h, "h", positive = TRUE)
  check_count(reps, "reps", minimum = 1)
  horizon <- length(risk)
  p_true <- true_odds_ratio * risk / (1 - risk + true_odds_ratio * risk)
  y <- matrix(as.numeric(stats::runif(reps * horizon) < rep(p_true, each = reps)), reps, horizon)
  W <- y * log(odds_ratio) - rep(log(1 - risk + odds_ratio * risk), each = reps)
  llr_paths(W, h)
}

#' Decision interval from a simulation
#'
#' The decision interval `h` such that a fraction `alpha` of the simulated **in-control** series have a CUSUM maximum above it: the
#' `1 - alpha` quantile of `max_stat` from [sim_cusum_poisson()] or [sim_cusum_risk()]. The false alarm probability over the length of
#' the series is then `alpha`.
#'
#' @param sim The result of [sim_cusum_poisson()] or [sim_cusum_risk()] with the true change set to 1.
#' @param alpha False alarm probability over the horizon of the series. Default 0.05.
#' @return A single number.
#' @export
cusum_limit <- function(sim, alpha = 0.05) {
  if (!is.data.frame(sim) || !"max_stat" %in% names(sim)) stop("`sim` must come from sim_cusum_poisson() or sim_cusum_risk().", call. = FALSE)
  check_probability(alpha, "alpha")
  unname(stats::quantile(sim$max_stat, 1 - alpha, type = 1))
}
