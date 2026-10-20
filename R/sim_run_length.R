# Draw `n_runs` runs of z-scores for the control levels.
draw_z <- function(n_runs, shift, dist, sdlog, scale = rep(1, length(shift))) {
  n_levels <- length(shift)
  if (dist == "normal") {
    z <- matrix(stats::rnorm(n_runs * n_levels), n_runs, n_levels)
  } else {
    # standardize the lognormal so that its true mean is 0 and its true SD is 1
    mu <- exp(sdlog^2 / 2)
    sigma <- sqrt((exp(sdlog^2) - 1) * exp(sdlog^2))
    z <- (matrix(stats::rlnorm(n_runs * n_levels, sdlog = sdlog), n_runs, n_levels) - mu) / sigma
  }
  z * rep(scale, each = n_runs) + rep(shift, each = n_runs)
}

# Run length (in runs) of one simulated sequence; NA if no rule fires within `max_runs` runs.
one_run_length <- function(rules, shift, max_runs, dist, sdlog, scale = rep(1, length(shift))) {
  have <- 0
  target <- min(max_runs, 500)
  z <- NULL
  repeat {
    z <- rbind(z, draw_z(target - have, shift, dist, sdlog, scale))
    fired <- rowSums(rule_fires(z, rules)) > 0
    if (any(fired)) return(which(fired)[1])
    if (target >= max_runs) return(NA_integer_)
    have <- target
    target <- min(max_runs, target * 4)
  }
}

#' Simulate the run length of a control rule
#'
#' Simulates control results run by run, applies the rules after every run, and records the
#' number of runs until the first rejection (the *run length*). The average of many run lengths
#' estimates the **average run length** (ARL). Use `shift = 0` for the in-control ARL and a
#' positive shift for the out-of-control ARL.
#'
#' Each run has one control result per level. The results are independent and, by default,
#' standard normal after standardizing (z-scores), shifted by `shift` standard deviations. With
#' `dist = "lognormal"` the results come from a lognormal distribution that is standardized to
#' have mean 0 and SD 1 *before* the shift, which shows what happens to the rules when the data
#' are skewed but the limits are still set at plus or minus 1, 2, and 3 SD of the true
#' distribution.
#'
#' The results can also have a standard deviation other than 1, set by `scale`. Together,
#' `shift` and `scale` describe a chart whose limits were **estimated**: if the true mean and SD
#' are \eqn{\mu}{mu} and \eqn{\sigma}{sigma} and the limits are placed at the estimates
#' \eqn{\hat\mu}{mu-hat} and \eqn{\hat\sigma}{sigma-hat}, the standardized results are
#' \eqn{(X - \hat\mu)/\hat\sigma}{(X - mu-hat)/sigma-hat}, which have mean shift
#' \eqn{(\mu - \hat\mu)/\hat\sigma}{(mu - mu-hat)/sigma-hat} and scale
#' \eqn{\sigma/\hat\sigma}{sigma/sigma-hat}.
#'
#' For a rule based on single independent points, such as `"1-3s"`, the run length is geometric
#' and \eqn{\text{ARL} = 1/p}{ARL = 1/p}. For rules that look at patterns (such as `"2-2s"`,
#' `"4-1s"`, and `"10x"`) and for combinations of rules, consecutive checks overlap, the run
#' length is not geometric, and simulation is the practical way to get the ARL.
#'
#' @param rules Rule names (see [apply_rules()]): any of `"1-2s"`, `"1-3s"`, `"2-2s"`, `"R-4s"`,
#'   `"4-1s"`, `"10x"`, or `"multirule"`. The rules are applied together: a run is rejected when any
#'   of them fires. The default is `"multirule"`.
#' @param shift Shift of the process mean in SD units. A single value is used for every level, or
#'   give one value per level.
#' @param levels Number of control levels per run. Default 2. Rules that compare levels within
#'   a run (R-4s, and the within-run form of 2-2s) need at least 2.
#' @param reps Number of simulated sequences. Default 2000.
#' @param max_runs Maximum length of a simulated sequence. Sequences with no rejection in
#'   `max_runs` runs are counted as `max_runs` and reported as censored. Default 5000.
#' @param dist Distribution of the controls: `"normal"` (default) or `"lognormal"`.
#' @param sdlog For `dist = "lognormal"`, the standard deviation of the logarithm. Default 0.5.
#' @param scale Standard deviation of the standardized results, relative to the SD used to set the
#'   limits. A single value is used for every level, or give one value per level. Default 1.
#'   A value above 1 means the process varies more than the limits assume.
#' @return An object of class `sqc_run_lengths`: a vector of `reps` run lengths, with a print
#'   method that shows the ARL, the standard deviation, and the quartiles. Use [arl()] for the ARL.
#' @examples
#' set.seed(1)
#' # In-control ARL of the 1-3s rule with two control levels per run:
#' # exact value 1 / (1 - (1 - 0.0027)^2) = 185
#' sim_run_length("1-3s", reps = 500)
#'
#' # A 2 SD shift is detected much faster by the multirule
#' sim_run_length("multirule", shift = 2, reps = 500)
#'
#' # Skewed (lognormal) data: the in-control ARL of one level is far below 370
#' sim_run_length("1-3s", levels = 1, dist = "lognormal", reps = 500)
#'
#' # Limits that are too narrow: the process SD is 1.5 times the SD used for the limits
#' # (exact ARL of 1-3s: 1 / (2 * (1 - pnorm(3 / 1.5))) = 22 runs)
#' sim_run_length("1-3s", levels = 1, scale = 1.5, reps = 500)
#' @seealso [arl()], [arl_table()], [apply_rules()]
#' @export
sim_run_length <- function(rules = "multirule", shift = 0, levels = 2, reps = 2000, max_runs = 5000,
                           dist = c("normal", "lognormal"), sdlog = 0.5, scale = 1) {
  rules <- expand_rules(rules)
  dist <- match.arg(dist)
  check_whole(levels, "levels", 1)
  check_whole(reps, "reps", 1)
  check_whole(max_runs, "max_runs", 2)
  if (!is.numeric(shift) || anyNA(shift) || !(length(shift) %in% c(1, levels))) {
    stop("`shift` must be a number, or one number per level.", call. = FALSE)
  }
  if (!is.numeric(sdlog) || length(sdlog) != 1 || sdlog <= 0) {
    stop("`sdlog` must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(scale) || anyNA(scale) || any(scale <= 0) || !(length(scale) %in% c(1, levels))) {
    stop("`scale` must be a positive number, or one positive number per level.", call. = FALSE)
  }
  if (levels < 2 && all(rules %in% "R-4s")) {
    stop("The R-4s rule needs at least 2 levels.", call. = FALSE)
  }
  shift <- rep_len(shift, levels)
  scale <- rep_len(scale, levels)
  rl <- vapply(seq_len(reps), function(i) one_run_length(rules, shift, max_runs, dist, sdlog, scale), numeric(1))
  censored <- sum(is.na(rl))
  rl[is.na(rl)] <- max_runs
  structure(as.integer(rl), class = "sqc_run_lengths",
            rules = rules, shift = shift, scale = scale, levels = levels, dist = dist,
            max_runs = max_runs, censored = censored)
}

check_whole <- function(x, name, minimum) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || x < minimum || x != round(x)) {
    stop("`", name, "` must be a whole number of at least ", minimum, ".", call. = FALSE)
  }
  invisible(x)
}

#' Average run length
#'
#' The average of the run lengths from [sim_run_length()].
#'
#' @param x An object returned by [sim_run_length()].
#' @return The average run length (a number).
#' @examples
#' set.seed(1)
#' arl(sim_run_length("1-2s", reps = 500))
#' @seealso [sim_run_length()], [arl_table()]
#' @export
arl <- function(x) {
  if (!inherits(x, "sqc_run_lengths")) stop("`x` must come from sim_run_length().", call. = FALSE)
  mean(unclass(x))
}

#' @export
print.sqc_run_lengths <- function(x, digits = 1, ...) {
  v <- unclass(x)
  q <- stats::quantile(v, c(0.25, 0.5, 0.75), names = FALSE)
  cat("<sqc_run_lengths> ", length(v), " simulated sequences\n", sep = "")
  cat("  rules:    ", paste(attr(x, "rules"), collapse = ", "), "\n", sep = "")
  cat("  shift:    ", paste(attr(x, "shift"), collapse = ", "), " SD (", attr(x, "levels"),
      " control level(s), ", attr(x, "dist"), ")\n", sep = "")
  sc <- attr(x, "scale")
  if (!is.null(sc) && any(sc != 1)) cat("  scale:    ", paste(sc, collapse = ", "), "\n", sep = "")
  cat("  ARL:      ", format(round(mean(v), digits), nsmall = digits), " runs (SD ",
      format(round(stats::sd(v), digits), nsmall = digits), ")\n", sep = "")
  cat("  quartiles:", paste(format(q), collapse = ", "), "\n", sep = " ")
  if (attr(x, "censored") > 0) {
    cat("  censored: ", attr(x, "censored"), " sequences had no rejection in ", attr(x, "max_runs"),
        " runs\n", sep = "")
  }
  invisible(x)
}

#' Table of average run lengths
#'
#' Estimates the average run length of several rules at several shifts with
#' [sim_run_length()] and returns the results as a tibble.
#'
#' @param rules A list of rule sets, one per row of the table. Each element is a character vector
#'   of rule names (see [apply_rules()]), and the names of the list label the rows. Unnamed elements
#'   are labeled with their rules joined by `+`. A plain character vector is treated as one
#'   rule set per element.
#' @param shifts Shifts of the process mean in SD units. Default `c(0, 1, 2)`.
#' @param ... Further arguments passed to [sim_run_length()] (`levels`, `reps`, `max_runs`,
#'   `dist`, `sdlog`, `scale`).
#' @return A tibble with columns `rule`, `shift`, `arl`, `sd` (the standard deviation of the run
#'   length), and `se` (the standard error of the ARL estimate).
#' @examples
#' set.seed(1)
#' arl_table(list("1-3s", "1-2s", Multirule = "multirule"), shifts = c(0, 2), reps = 300)
#' @seealso [sim_run_length()]
#' @export
arl_table <- function(rules, shifts = c(0, 1, 2), ...) {
  if (is.character(rules)) rules <- as.list(rules)
  labels <- names(rules)
  if (is.null(labels)) labels <- rep("", length(rules))
  labels <- ifelse(labels == "", vapply(rules, paste, character(1), collapse = "+"), labels)
  out <- list()
  for (i in seq_along(rules)) {
    for (s in shifts) {
      rl <- sim_run_length(rules[[i]], shift = s, ...)
      v <- unclass(rl)
      out[[length(out) + 1]] <- tibble::tibble(rule = labels[i], shift = s, arl = mean(v),
                                               sd = stats::sd(v), se = stats::sd(v) / sqrt(length(v)))
    }
  }
  do.call(rbind, out)
}
