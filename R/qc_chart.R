#' Compute a Shewhart control chart
#'
#' Computes the plotted statistics, center line, control limits, and signals for an
#' \eqn{\bar X}{Xbar} or \eqn{R} chart from subgroup data. Draw the chart with [autoplot.qc_chart()].
#'
#' The center line and limits are estimated from the data unless `center` or `sigma` are given.
#' With subgroup size \eqn{n}, the within-subgroup standard deviation is estimated by
#' \eqn{\hat\sigma = \bar R / d_2}{sigma-hat = Rbar / d_2}, and
#' * the \eqn{\bar X}{Xbar} chart has center line \eqn{\bar{\bar x}}{xbarbar} (or `center`) and limits
#'   \eqn{\text{CL} \pm k\hat\sigma/\sqrt{n}}{CL +/- k*sigma-hat/sqrt(n)};
#' * the \eqn{R} chart has center line \eqn{\bar R = d_2\hat\sigma}{Rbar = d_2*sigma-hat} and limits
#'   \eqn{\max(0, \bar R - k d_3\hat\sigma)}{max(0, Rbar - k d_3*sigma-hat)} and \eqn{\bar R + k d_3\hat\sigma}{Rbar + k d_3*sigma-hat}.
#'
#' Two signals are flagged: a point beyond the limits, and a run of `run_length` or more
#' consecutive points on one side of the center line (the points from the `run_length`th
#' onward are flagged, since that is when the run is recognized).
#'
#' @param data Subgroup data, either a matrix with one row per subgroup (sample) and one column
#'   per observation, or a data frame with one row per observation (use `value` and `sample`).
#'   All subgroups must have the same size \eqn{n \ge 2}{n >= 2}. Subgroups are plotted in the order
#'   in which they first appear.
#' @param type `"xbar"` for the \eqn{\bar X}{Xbar} chart or `"R"` for the \eqn{R} chart (case is ignored).
#' @param value,sample For a data frame: the column holding the measurements and the column
#'   identifying the subgroup (bare names, as in dplyr).
#' @param center,sigma Optional known (or target) center line and within-subgroup standard
#'   deviation \eqn{\sigma}{sigma}. Supplying them gives limits that are not re-estimated, which is
#'   how a Phase II chart is built from Phase I estimates.
#' @param k Number of standard deviations from the center line to the limits: a positive
#'   number. Default 3.
#' @param run_length Length of a one-sided run that counts as a signal: a whole number of at
#'   least 2, or `Inf` to turn the run signal off. Default 7.
#' @section Errors and warnings:
#' Invalid arguments give an error: `k` must be a positive number, `sigma` a positive finite
#' number, `center` a finite number (positive for an \eqn{R} chart), and `run_length` a whole
#' number of at least 2 (or `Inf`). Subgroups must all have the same size, at least 2, with no
#' missing values. A warning is given if the estimated \eqn{\sigma}{sigma} is 0 (every subgroup
#' has range 0), because the limits then collapse onto the center line.
#'
#' @return An object of class `qc_chart`, a list with elements `type`, `statistics`, `center`,
#'   `limits` (`lcl`, `ucl`), `sigma`, `n` (subgroup size), `m` (number of subgroups), `k`,
#'   `run_length`, and `signals` (indices of points beyond the limits, `beyond`, and in a
#'   run, `run`).
#' @examples
#' # A matrix with one row per sample
#' set.seed(1)
#' x <- matrix(rnorm(20 * 4, mean = 500, sd = 1), ncol = 4)
#' q <- qc_chart(x, type = "xbar")
#' q
#'
#' # A data frame with one row per bottle
#' qc_chart(filler, type = "r", value = fill_ml, sample = cycle)
#'
#' # Phase II: monitor new data with limits fixed from Phase I
#' phase1 <- qc_chart(x, type = "xbar")
#' new <- matrix(rnorm(5 * 4, mean = 501, sd = 1), ncol = 4)
#' qc_chart(new, type = "xbar", center = phase1$center, sigma = phase1$sigma)
#' @export
qc_chart <- function(data, type = c("xbar", "R"), value = NULL, sample = NULL,
                     center = NULL, sigma = NULL, k = 3, run_length = 7) {
  type <- match.arg(tolower(type[1]), c("xbar", "r"))
  value <- rlang::enquo(value)
  sample <- rlang::enquo(sample)
  check_number(k, "k", positive = TRUE)
  check_run_length(run_length)
  if (!is.null(sigma)) check_number(sigma, "sigma", positive = TRUE)
  if (!is.null(center)) check_number(center, "center", positive = (type == "r"))
  x <- subgroup_matrix(data, value, sample)
  n <- ncol(x)
  m <- nrow(x)

  columns <- asplit(x, 2)
  ranges <- as.vector(do.call(pmax, columns) - do.call(pmin, columns))   # row-wise max - min, vectorized
  d2n <- d2(n)
  if (type == "xbar") {
    statistics <- rowMeans(x)
    sigma_hat <- if (is.null(sigma)) mean(ranges) / d2n else sigma
    center_line <- if (is.null(center)) mean(statistics) else center
    half <- k * sigma_hat / sqrt(n)
    limits <- c(lcl = center_line - half, ucl = center_line + half)
  } else {
    statistics <- ranges
    sigma_hat <- if (!is.null(sigma)) sigma else if (!is.null(center)) center / d2n else mean(ranges) / d2n
    center_line <- if (!is.null(center)) center else d2n * sigma_hat
    half <- k * d3(n) * sigma_hat
    limits <- c(lcl = max(0, center_line - half), ucl = center_line + half)
  }

  if (sigma_hat == 0) {
    warning("The estimated sigma is 0 (every subgroup has range 0), so the limits collapse onto ",
            "the center line.", call. = FALSE)
  }

  side <- ifelse(statistics > center_line, 1, -1)
  runs <- rle(side)
  in_run <- sequence(runs$lengths) >= run_length   # the signal fires at the `run_length`th point of a run
  structure(
    list(type = type, statistics = statistics, center = center_line, limits = limits,
         sigma = sigma_hat, n = n, m = m, k = k, run_length = run_length,
         signals = list(beyond = which(statistics < limits[["lcl"]] | statistics > limits[["ucl"]]),
                        run = which(in_run))),
    class = "qc_chart")
}

# `x` must be a single finite number (and positive when `positive = TRUE`).
check_number <- function(x, name, positive = FALSE) {
  ok <- is.numeric(x) && length(x) == 1 && is.finite(x) && (!positive || x > 0)
  if (!ok) {
    stop("`", name, "` must be a single ", if (positive) "positive ", "finite number.", call. = FALSE)
  }
  invisible(x)
}

# `run_length` must be a whole number of at least 2, or Inf (no run signal).
check_run_length <- function(x) {
  ok <- is.numeric(x) && length(x) == 1 && !is.na(x) && x >= 2 && (is.infinite(x) || x == round(x))
  if (!ok) {
    stop("`run_length` must be a whole number of at least 2 (or `Inf` to turn the run signal off).",
         call. = FALSE)
  }
  invisible(x)
}

# Turn matrix or data-frame input into an m x n matrix (one row per subgroup).
subgroup_matrix <- function(data, value, sample) {
  if (is.data.frame(data)) {
    if (rlang::quo_is_null(value) || rlang::quo_is_null(sample)) {
      stop("For a data frame, supply both `value` and `sample` (the measurement and the subgroup columns). ",
           "If each row is already a subgroup (wide format), convert it with `as.matrix()`.",
           call. = FALSE)
    }
    v <- rlang::eval_tidy(value, data)
    s <- rlang::eval_tidy(sample, data)
    groups <- split(v, factor(s, levels = unique(s)))
  } else if (is.matrix(data)) {
    if (!rlang::quo_is_null(value) || !rlang::quo_is_null(sample)) {
      stop("`value` and `sample` are only used with a data frame.", call. = FALSE)
    }
    groups <- split(data, row(data))
  } else {
    stop("`data` must be a matrix (one row per subgroup) or a data frame.", call. = FALSE)
  }
  sizes <- lengths(groups)
  if (length(unique(sizes)) != 1) {
    stop("All subgroups must have the same size (found sizes ", paste(unique(sizes), collapse = ", "), ").",
         call. = FALSE)
  }
  if (sizes[[1]] < 2) {
    stop("Subgroups must have at least 2 observations.", call. = FALSE)
  }
  x <- do.call(rbind, unname(groups))
  if (!is.numeric(x) || anyNA(x)) {
    stop("Measurements must be numeric with no missing values.", call. = FALSE)
  }
  x
}

#' @export
print.qc_chart <- function(x, digits = 3, ...) {
  name <- if (x$type == "xbar") "X-bar chart" else "R chart"
  cat("<qc_chart> ", name, ": ", x$m, " subgroups of size ", x$n, "\n", sep = "")
  cat("  center line:      ", format(round(x$center, digits), nsmall = digits), "\n", sep = "")
  cat("  limits (LCL, UCL): ", format(round(x$limits[["lcl"]], digits), nsmall = digits), ", ",
      format(round(x$limits[["ucl"]], digits), nsmall = digits), "\n", sep = "")
  cat("  sigma-hat:        ", format(round(x$sigma, digits), nsmall = digits), "\n", sep = "")
  cat("  signals:          ", length(x$signals$beyond), " beyond the limits; ",
      length(x$signals$run), " in runs of ", x$run_length, "+\n", sep = "")
  invisible(x)
}

#' @export
as.data.frame.qc_chart <- function(x, ...) {
  status <- rep("No signal", x$m)
  status[x$signals$run] <- "Run"
  status[x$signals$beyond] <- "Beyond limits"
  data.frame(sample = seq_len(x$m), statistic = x$statistics, status = status)
}
