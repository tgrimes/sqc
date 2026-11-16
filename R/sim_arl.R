#' Simulate the run length of a Shewhart, CUSUM, or EWMA chart
#'
#' Simulates independent normal results with a shift of the mean of `shift` standard deviations and records the number of
#' results until the chart signals (the *run length*). The average of many run lengths estimates the ARL, which [shewhart_arl()],
#' [cusum_arl()], and [ewma_arl()] compute numerically; the simulation also gives the distribution of the run length.
#'
#' The shift starts at the first result. The Shewhart chart has limits at \eqn{\pm L}; the CUSUM is the tabular CUSUM of
#' [cusum_chart()] (`sides` 1 or 2, started at `head_start`); the EWMA uses the asymptotic limits at
#' \eqn{\pm L\sqrt{\lambda/(2-\lambda)}}{+/- L sqrt(lambda/(2 - lambda))}.
#'
#' @param chart `"shewhart"`, `"cusum"`, or `"ewma"`.
#' @param shift Shift of the mean in standard deviations. Default 0.
#' @param reps Number of simulated sequences. Default 2000.
#' @param max_runs Longest sequence simulated; sequences without a signal count as `max_runs` and are reported as censored. Default 20000.
#' @param L Width of the limits (Shewhart and EWMA). Default 3.
#' @param k,h,sides,head_start CUSUM parameters (see [cusum_chart()]).
#' @param lambda EWMA weight. Default 0.2.
#' @return An object of class `sqc_chart_run_lengths`: a vector of `reps` run lengths with a print method, to be used with [arl()].
#' @seealso [shewhart_arl()], [cusum_arl()], [ewma_arl()], [sim_run_length()]
#' @examples
#' set.seed(1)
#' sim_arl("cusum", shift = 1, k = 0.5, h = 5, reps = 1000)
#' cusum_arl(1, k = 0.5, h = 5)     # numerical value
#' @export
sim_arl <- function(chart = c("shewhart", "cusum", "ewma"), shift = 0, reps = 2000, max_runs = 20000, L = 3,
                    k = 0.5, h = 5, sides = 2, head_start = 0, lambda = 0.2) {
  chart <- match.arg(chart)
  check_finite(shift, "shift"); check_count(reps, "reps", minimum = 1); check_count(max_runs, "max_runs", minimum = 2)
  check_number(L, "L", positive = TRUE)
  if (chart == "cusum") {
    check_number(k, "k", positive = TRUE); check_number(h, "h", positive = TRUE)
    if (!sides %in% c(1, 2)) stop("`sides` must be 1 or 2.", call. = FALSE)
    if (head_start < 0 || head_start >= h) stop("`head_start` must be between 0 and `h`.", call. = FALSE)
  }
  if (chart == "ewma" && !(is.numeric(lambda) && length(lambda) == 1 && lambda > 0 && lambda <= 1)) {
    stop("`lambda` must be a single number in (0, 1].", call. = FALSE)
  }
  rl <- rep(NA_integer_, reps)
  active <- seq_len(reps)
  up <- dn <- rep(head_start, reps)
  zz <- numeric(reps)
  cc <- L * sqrt(lambda / (2 - lambda))
  for (t in seq_len(max_runs)) {
    if (!length(active)) break
    x <- stats::rnorm(length(active), shift)
    hit <- switch(chart,
      shewhart = abs(x) > L,
      cusum = {
        up[active] <- pmax(0, up[active] + x - k)
        dn[active] <- pmax(0, dn[active] - x - k)
        if (sides == 2) (up[active] > h) | (dn[active] > h) else up[active] > h
      },
      ewma = {
        zz[active] <- lambda * x + (1 - lambda) * zz[active]
        abs(zz[active]) > cc
      })
    rl[active[hit]] <- t
    active <- active[!hit]
  }
  censored <- sum(is.na(rl))
  rl[is.na(rl)] <- max_runs
  structure(as.integer(rl), class = "sqc_chart_run_lengths", chart = chart, shift = shift, max_runs = max_runs, censored = censored,
            params = switch(chart, shewhart = list(L = L), cusum = list(k = k, h = h, sides = sides, head_start = head_start),
                            ewma = list(lambda = lambda, L = L)))
}

#' @export
print.sqc_chart_run_lengths <- function(x, digits = 1, ...) {
  v <- unclass(x)
  q <- stats::quantile(v, c(0.25, 0.5, 0.75), names = FALSE)
  p <- attr(x, "params")
  cat("<sqc_chart_run_lengths> ", length(v), " simulated sequences\n", sep = "")
  cat("  chart:    ", attr(x, "chart"), " (", paste(names(p), unlist(p), sep = " = ", collapse = ", "), ")\n", sep = "")
  cat("  shift:    ", attr(x, "shift"), " SD\n", sep = "")
  cat("  ARL:      ", format(round(mean(v), digits), nsmall = digits), " (SD ", format(round(stats::sd(v), digits), nsmall = digits), ")\n", sep = "")
  cat("  quartiles:", paste(format(q), collapse = ", "), "\n", sep = " ")
  if (attr(x, "censored") > 0) cat("  censored: ", attr(x, "censored"), " sequences had no signal in ", attr(x, "max_runs"), " results\n", sep = "")
  invisible(x)
}
