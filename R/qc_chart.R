#' Compute a Shewhart control chart
#'
#' Computes the plotted statistics, center line, control limits, and signals for an
#' \eqn{\bar X}{Xbar} or \eqn{R} chart from subgroup data, or for an individuals (I) or moving range
#' (MR) chart from individual results. Draw the chart with [autoplot.qc_chart()].
#'
#' The center line and limits are estimated from the data unless `center` or `sigma` are given.
#' With subgroup size \eqn{n}, the within-subgroup standard deviation is estimated by
#' \eqn{\hat\sigma = \bar R / d_2}{sigma-hat = Rbar / d_2}, and
#' * the \eqn{\bar X}{Xbar} chart has center line \eqn{\bar{\bar x}}{xbarbar} (or `center`) and limits
#'   \eqn{\text{CL} \pm k\hat\sigma/\sqrt{n}}{CL +/- k*sigma-hat/sqrt(n)};
#' * the \eqn{R} chart has center line \eqn{\bar R = d_2\hat\sigma}{Rbar = d_2*sigma-hat} and limits
#'   \eqn{\max(0, \bar R - k d_3\hat\sigma)}{max(0, Rbar - k d_3*sigma-hat)} and \eqn{\bar R + k d_3\hat\sigma}{Rbar + k d_3*sigma-hat}.
#'
#' For **individual results** (one observation per time point) the "subgroups" are consecutive
#' pairs. The moving range of the \eqn{i}th result is \eqn{MR_i = |x_i - x_{i-1}|}{MR_i = |x_i - x_(i-1)|},
#' \eqn{\hat\sigma = \overline{MR} / d_2}{sigma-hat = MRbar / d_2} with \eqn{d_2 = 1.128} (the value for
#' \eqn{n = 2}), and
#' * the I chart plots the results with center line \eqn{\bar x}{xbar} (or `center`) and limits
#'   \eqn{\text{CL} \pm k\hat\sigma}{CL +/- k*sigma-hat} (for \eqn{k = 3}, \eqn{\bar x \pm 2.66\,\overline{MR}}{xbar +/- 2.66 MRbar});
#' * the MR chart plots the moving ranges (starting at the second result) with the R chart limits for
#'   \eqn{n = 2} (for \eqn{k = 3}, \eqn{0} and \eqn{3.267\,\overline{MR}}{3.267 MRbar}).
#'
#' A point beyond the limits is always flagged as a signal. Patterns within the limits are not
#' checked unless you ask for them: set `run_length` to flag a run of that many (or more)
#' consecutive points on one side of the center line. The points from the `run_length`th onward
#' are flagged, since that is when the run is recognized.
#'
#' In Phase I, subgroups or results with an assignable cause are left out of the estimation with
#' `exclude`. They are still plotted (and flagged if they lie beyond the limits), but they do not
#' contribute to the center line or to \eqn{\hat\sigma}{sigma-hat}. For individual results a moving range is
#' left out when either of its two results is excluded.
#'
#' @param data Subgroup data, either a matrix with one row per subgroup (sample) and one column
#'   per observation, or a data frame with one row per observation (use `value` and `sample`).
#'   All subgroups must have the same size \eqn{n \ge 2}{n >= 2}. Subgroups are plotted in the order
#'   in which they first appear. For `type = "I"` or `"MR"`, a numeric vector of individual results
#'   in time order, or a data frame with the results in the column `value`.
#' @param type `"xbar"` for the \eqn{\bar X}{Xbar} chart, `"R"` for the \eqn{R} chart, `"I"` for the
#'   individuals chart, or `"MR"` for the moving range chart (case is ignored).
#' @param value,sample For a data frame: the column holding the measurements and the column
#'   identifying the subgroup (bare names, as in dplyr). The I and MR charts use `value` only.
#' @param center,sigma Optional known (or target) center line and within-subgroup standard
#'   deviation \eqn{\sigma}{sigma}. Supplying them gives limits that are not re-estimated, which is
#'   how a Phase II chart is built from Phase I estimates.
#' @param k Number of standard deviations from the center line to the limits: a positive
#'   number. Default 3.
#' @param run_length Length of a one-sided run that counts as a signal: a whole number of at
#'   least 2. The default, `NULL`, does not look for runs (`Inf` has the same effect).
#' @param exclude Positions (subgroup numbers, or the observation numbers for I and MR charts) to
#'   leave out of the estimation of the center line and \eqn{\sigma}{sigma}, for example runs with a known
#'   assignable cause. Default `NULL`.
#' @section Errors and warnings:
#' Invalid arguments give an error: `k` must be a positive number, `sigma` a positive finite
#' number, `center` a finite number (positive for an \eqn{R} or MR chart), `run_length` `NULL`, `Inf`,
#' or a whole number of at least 2, and `exclude` whole numbers between 1 and the number of
#' subgroups (or results), leaving enough to estimate from. Subgroups must all have the same size, at
#' least 2, with no missing values. A warning is given if the estimated \eqn{\sigma}{sigma} is 0
#' (every subgroup has range 0), because the limits then collapse onto the center line.
#'
#' @return An object of class `qc_chart`, a list with elements `type`, `statistics`, `index` (the
#'   observation or subgroup number of each plotted statistic: `1:m`, or `2:m` for the MR chart),
#'   `center`, `limits` (`lcl`, `ucl`), `sigma`, `n` (subgroup size, 1 for the I chart and 2 for the
#'   MR chart), `m` (number of subgroups, or of individual results), `k`, `run_length` (`Inf` when
#'   runs are not checked), `exclude`, and `signals` (positions of points beyond the limits,
#'   `beyond`, and in a run, `run`, as subgroup or observation numbers; the latter is empty unless
#'   `run_length` is set).
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
#' # Also flag runs of 7 or more points on one side of the center line (off by default)
#' qc_chart(x, type = "xbar", run_length = 7)
#'
#' # Phase II: monitor new data with limits fixed from Phase I
#' phase1 <- qc_chart(x, type = "xbar")
#' new <- matrix(rnorm(5 * 4, mean = 501, sd = 1), ncol = 4)
#' qc_chart(new, type = "xbar", center = phase1$center, sigma = phase1$sigma)
#'
#' # Individual results: the I and MR charts, with sigma estimated from the moving ranges
#' y <- c(rnorm(15, mean = 100, sd = 3), 118, rnorm(14, mean = 100, sd = 3))
#' qc_chart(y, type = "I")
#' qc_chart(y, type = "MR")
#'
#' # Leave a run with a known assignable cause (run 16) out of the estimation
#' qc_chart(y, type = "I", exclude = 16)
#' @export
qc_chart <- function(data, type = c("xbar", "R", "I", "MR"), value = NULL, sample = NULL,
                     center = NULL, sigma = NULL, k = 3, run_length = NULL, exclude = NULL) {
  type <- match.arg(tolower(type[1]), c("xbar", "r", "i", "mr"))
  value <- rlang::enquo(value)
  sample <- rlang::enquo(sample)
  check_number(k, "k", positive = TRUE)
  check_run_length(run_length)
  if (is.null(run_length)) run_length <- Inf   # no run signal
  if (!is.null(sigma)) check_number(sigma, "sigma", positive = TRUE)
  if (!is.null(center)) check_number(center, "center", positive = (type %in% c("r", "mr")))

  if (type %in% c("i", "mr")) {
    x <- individual_values(data, value, sample)
    m <- length(x)
    use <- use_index(exclude, m, minimum = 3)
    mr <- c(NA, abs(diff(x)))
    mr_in_use <- !is.na(mr) & use & c(FALSE, use[-m])      # both results of the pair are in use
    d2n <- d2(2)
    if (type == "i") {
      statistics <- x
      index <- seq_len(m)
      sigma_hat <- if (is.null(sigma)) mean(mr[mr_in_use]) / d2n else sigma
      center_line <- if (is.null(center)) mean(x[use]) else center
      half <- k * sigma_hat
      limits <- c(lcl = center_line - half, ucl = center_line + half)
      n <- 1
    } else {
      statistics <- mr[-1]
      index <- 2:m
      sigma_hat <- if (!is.null(sigma)) sigma else if (!is.null(center)) center / d2n else mean(mr[mr_in_use]) / d2n
      center_line <- if (!is.null(center)) center else d2n * sigma_hat
      half <- k * d3(2) * sigma_hat
      limits <- c(lcl = max(0, center_line - half), ucl = center_line + half)
      n <- 2
    }
    warn_zero <- "The estimated sigma is 0 (all moving ranges are 0), so the limits collapse onto the center line."
  } else {
    x <- subgroup_matrix(data, value, sample)
    n <- ncol(x)
    m <- nrow(x)
    use <- use_index(exclude, m, minimum = 1)
    columns <- asplit(x, 2)
    ranges <- as.vector(do.call(pmax, columns) - do.call(pmin, columns))   # row-wise max - min, vectorized
    d2n <- d2(n)
    index <- seq_len(m)
    if (type == "xbar") {
      statistics <- rowMeans(x)
      sigma_hat <- if (is.null(sigma)) mean(ranges[use]) / d2n else sigma
      center_line <- if (is.null(center)) mean(statistics[use]) else center
      half <- k * sigma_hat / sqrt(n)
      limits <- c(lcl = center_line - half, ucl = center_line + half)
    } else {
      statistics <- ranges
      sigma_hat <- if (!is.null(sigma)) sigma else if (!is.null(center)) center / d2n else mean(ranges[use]) / d2n
      center_line <- if (!is.null(center)) center else d2n * sigma_hat
      half <- k * d3(n) * sigma_hat
      limits <- c(lcl = max(0, center_line - half), ucl = center_line + half)
    }
    warn_zero <- "The estimated sigma is 0 (every subgroup has range 0), so the limits collapse onto the center line."
  }

  if (sigma_hat == 0) warning(warn_zero, call. = FALSE)

  side <- ifelse(statistics > center_line, 1, -1)
  runs <- rle(side)
  in_run <- sequence(runs$lengths) >= run_length   # the signal fires at the `run_length`th point of a run
  structure(
    list(type = type, statistics = statistics, index = index, center = center_line, limits = limits,
         sigma = sigma_hat, n = n, m = m, k = k, run_length = run_length,
         exclude = sort(unique(as.integer(exclude))),
         signals = list(beyond = index[which(statistics < limits[["lcl"]] | statistics > limits[["ucl"]])],
                        run = index[which(in_run)])),
    class = "qc_chart")
}

# Logical vector of length `m`: TRUE for the subgroups (or results) used in the estimation.
use_index <- function(exclude, m, minimum) {
  use <- rep(TRUE, m)
  if (is.null(exclude)) return(use)
  ok <- is.numeric(exclude) && !anyNA(exclude) && all(exclude == round(exclude)) &&
    all(exclude >= 1 & exclude <= m)
  if (!ok) {
    stop("`exclude` must be whole numbers between 1 and ", m, ".", call. = FALSE)
  }
  use[exclude] <- FALSE
  if (sum(use) < minimum) {
    stop("Too few subgroups (or results) are left to estimate from after `exclude`.", call. = FALSE)
  }
  use
}

# Individual results as a numeric vector (from a vector, or from the `value` column of a data frame).
individual_values <- function(data, value, sample) {
  if (is.data.frame(data)) {
    if (rlang::quo_is_null(value)) {
      stop("For a data frame, supply `value` (the column of individual results).", call. = FALSE)
    }
    x <- rlang::eval_tidy(value, data)
  } else if (is.numeric(data) && is.null(dim(data))) {
    if (!rlang::quo_is_null(value)) stop("`value` is only used with a data frame.", call. = FALSE)
    x <- data
  } else {
    stop("For an I or MR chart, `data` must be a numeric vector or a data frame with a `value` column.",
         call. = FALSE)
  }
  if (!rlang::quo_is_null(sample)) {
    stop("`sample` is not used for the I and MR charts.", call. = FALSE)
  }
  if (!is.numeric(x) || anyNA(x) || length(x) < 3) {
    stop("Individual results must be numeric with no missing values, at least 3 of them.", call. = FALSE)
  }
  x
}

# `x` must be a single finite number (and positive when `positive = TRUE`).
check_number <- function(x, name, positive = FALSE) {
  ok <- is.numeric(x) && length(x) == 1 && is.finite(x) && (!positive || x > 0)
  if (!ok) {
    stop("`", name, "` must be a single ", if (positive) "positive ", "finite number.", call. = FALSE)
  }
  invisible(x)
}

# `run_length` must be NULL or Inf (no run signal), or a whole number of at least 2.
check_run_length <- function(x) {
  ok <- is.null(x) ||
    (is.numeric(x) && length(x) == 1 && !is.na(x) && x >= 2 && (is.infinite(x) || x == round(x)))
  if (!ok) {
    stop("`run_length` must be a whole number of at least 2 (or `NULL` to not look for runs).",
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
  name <- c(xbar = "X-bar chart", r = "R chart", i = "I chart", mr = "MR chart")[[x$type]]
  what <- if (x$type %in% c("i", "mr")) paste0(x$m, " individual results") else paste0(x$m, " subgroups of size ", x$n)
  cat("<qc_chart> ", name, ": ", what, "\n", sep = "")
  cat("  center line:      ", format(round(x$center, digits), nsmall = digits), "\n", sep = "")
  cat("  limits (LCL, UCL): ", format(round(x$limits[["lcl"]], digits), nsmall = digits), ", ",
      format(round(x$limits[["ucl"]], digits), nsmall = digits), "\n", sep = "")
  cat("  sigma-hat:        ", format(round(x$sigma, digits), nsmall = digits), "\n", sep = "")
  if (length(x$exclude) > 0) {
    cat("  left out of the estimation: ", paste(x$exclude, collapse = ", "), "\n", sep = "")
  }
  cat("  signals:          ", length(x$signals$beyond), " beyond the limits", sep = "")
  if (is.finite(x$run_length)) {
    cat("; ", length(x$signals$run), " in runs of ", x$run_length, "+", sep = "")
  }
  cat("\n")
  invisible(x)
}

#' @export
as.data.frame.qc_chart <- function(x, ...) {
  status <- rep("No signal", length(x$statistics))
  status[match(x$signals$run, x$index)] <- "Run"
  status[match(x$signals$beyond, x$index)] <- "Beyond limits"
  data.frame(sample = x$index, statistic = x$statistics, status = status,
             excluded = x$index %in% x$exclude)
}
