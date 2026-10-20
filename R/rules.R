# The rule names. "1-2s" is a warning rule in Westgard's scheme; it is available here as a
# rejection rule so that its false-rejection rate can be studied.
rule_names <- c("1-2s", "1-3s", "2-2s", "R-4s", "4-1s", "10x")
multirule_names <- c("1-3s", "2-2s", "R-4s", "4-1s", "10x")

# Expand "multirule" and check that every rule name is known.
expand_rules <- function(rules) {
  if (!is.character(rules) || length(rules) == 0) {
    stop("`rules` must be a character vector of rule names.", call. = FALSE)
  }
  rules <- unlist(lapply(rules, function(r) if (identical(r, "multirule")) multirule_names else r))
  bad <- setdiff(rules, rule_names)
  if (length(bad) > 0) {
    stop("Unknown rule(s): ", paste(bad, collapse = ", "), ". The rules are ",
         paste(rule_names, collapse = ", "), ', or "multirule".', call. = FALSE)
  }
  unique(rules)
}

# TRUE at position t when the last k values of the logical vector `x` are all TRUE.
roll_all <- function(x, k) {
  n <- length(x)
  if (n < k) return(rep(FALSE, n))
  s <- as.numeric(stats::filter(as.numeric(x), rep(1, k), sides = 1))
  !is.na(s) & s == k
}

# TRUE where the condition holds for at least one level (columns of the logical matrix `m`).
any_level <- function(m, k) {
  Reduce(`|`, lapply(seq_len(ncol(m)), function(l) roll_all(m[, l], k)))
}

# Logical matrix (runs x rules): does each rule fire at each run, using the data up to that run?
# `z` is a matrix of z-scores with one row per run and one column per control level.
# Rules that compare levels within a run (R-4s, and the within-run form of 2-2s) need 2 or more levels.
rule_fires <- function(z, rules) {
  z <- as.matrix(z)
  m <- nrow(z)
  n_levels <- ncol(z)
  out <- matrix(FALSE, m, length(rules), dimnames = list(NULL, rules))
  above2 <- z > 2
  below2 <- z < -2
  for (r in rules) {
    out[, r] <- switch(r,
      "1-2s" = rowSums(abs(z) > 2) > 0,
      "1-3s" = rowSums(abs(z) > 3) > 0,
      "2-2s" = {
        within_run <- if (n_levels >= 2) rowSums(above2) >= 2 | rowSums(below2) >= 2 else rep(FALSE, m)
        within_run | any_level(above2, 2) | any_level(below2, 2)
      },
      "R-4s" = if (n_levels >= 2) rowSums(above2) >= 1 & rowSums(below2) >= 1 else rep(FALSE, m),
      "4-1s" = any_level(z > 1, 4) | any_level(z < -1, 4),
      "10x" = any_level(z > 0, 10) | any_level(z < 0, 10)
    )
  }
  out
}

#' Westgard rules: apply rejection rules to control z-scores
#'
#' Applies the Westgard rules to the standardized control results of a laboratory, with one row
#' per run and one column per control level. Each control result is first converted to a
#' z-score, \eqn{z = (x - \text{mean})/\text{SD}}{z = (x - mean)/SD}, using the established mean
#' and SD of that control material. A run is *rejected* when a rule fires.
#'
#' @section The rules:
#' Rules are written as \eqn{a_L}{a-Ls}: `a` control observations beyond `L` SD.
#'
#' | Rule | The rule fires when | Null probability per point |
#' |:--|:--|:--|
#' | `"1-2s"` | one control is beyond 2 SD (a warning rule in Westgard's scheme) | \eqn{2(1-\Phi(2)) = 0.0455}{2(1 - Phi(2)) = 0.0455} per control |
#' | `"1-3s"` | one control is beyond 3 SD | \eqn{2(1-\Phi(3)) = 0.0027}{2(1 - Phi(3)) = 0.0027} per control |
#' | `"2-2s"` | two consecutive controls are beyond the same 2 SD limit: two levels in one run, or one level in two consecutive runs | \eqn{2(1-\Phi(2))^2 = 0.00104}{2(1 - Phi(2))^2 = 0.00104} |
#' | `"R-4s"` | in one run, one control is above +2 SD and another is below -2 SD (needs 2 or more levels) | \eqn{2(1-\Phi(2))^2 = 0.00104}{2(1 - Phi(2))^2 = 0.00104} for two levels |
#' | `"4-1s"` | four consecutive controls of one level are beyond the same 1 SD limit | \eqn{2(1-\Phi(1))^4 = 0.00127}{2(1 - Phi(1))^4 = 0.00127} |
#' | `"10x"` | ten consecutive controls of one level are on the same side of the mean | \eqn{2(1/2)^{10} = 0.00195}{2(1/2)^10 = 0.00195} |
#'
#' `"multirule"` stands for `"1-3s"`, `"2-2s"`, `"R-4s"`, `"4-1s"`, and `"10x"` together. The null
#' probabilities assume independent standard normal controls (the process is in control).
#'
#' @section Context:
#' The Levey-Jennings chart and the Westgard rules are the standard vocabulary of clinical
#' laboratory quality control. In general SPC, the same ideas are known as an individuals chart
#' with sigma zones and the Western Electric and Nelson run rules.
#'
#' @param z A numeric matrix (or data frame) of z-scores with one row per run and one column per
#'   control level, or a numeric vector for a single level.
#' @param rules Rule names: any of `"1-2s"`, `"1-3s"`, `"2-2s"`, `"R-4s"`, `"4-1s"`, `"10x"`, or
#'   `"multirule"` (the default).
#' @return A tibble with one row per rule violation and the columns `run`, `rule`, and `level`
#'   (`NA` when the rule compares levels within a run).
#' @examples
#' z <- cbind(c(0.2, -0.5, 1.1, 2.3, 2.6, 0.1), c(0.4, 1.2, -0.3, 2.2, 0.5, -2.4))
#' apply_rules(z)
#' apply_rules(z, rules = "1-2s")
#' @seealso [lj_chart()] for a chart that applies the rules, and [sim_run_length()] for the
#'   average run length of a rule.
#' @export
apply_rules <- function(z, rules = "multirule") {
  rules <- expand_rules(rules)
  z <- as.matrix(z)
  if (!is.numeric(z) || anyNA(z)) stop("`z` must be numeric with no missing values.", call. = FALSE)
  m <- nrow(z)
  n_levels <- ncol(z)
  above2 <- z > 2
  below2 <- z < -2
  rows <- list()
  add <- function(run, rule, level) {
    if (length(run) > 0) rows[[length(rows) + 1]] <<- tibble::tibble(run = run, rule = rule, level = level)
  }
  for (r in rules) {
    if (r %in% c("1-2s", "1-3s")) {
      limit <- if (r == "1-2s") 2 else 3
      for (l in seq_len(n_levels)) {
        run <- which(abs(z[, l]) > limit)
        add(run, r, rep(l, length(run)))
      }
    } else if (r == "2-2s") {
      if (n_levels >= 2) {
        run <- which(rowSums(above2) >= 2 | rowSums(below2) >= 2)
        add(run, r, rep(NA_integer_, length(run)))
      }
      for (l in seq_len(n_levels)) {
        run <- which(roll_all(above2[, l], 2) | roll_all(below2[, l], 2))
        add(run, r, rep(l, length(run)))
      }
    } else if (r == "R-4s") {
      if (n_levels >= 2) {
        run <- which(rowSums(above2) >= 1 & rowSums(below2) >= 1)
        add(run, r, rep(NA_integer_, length(run)))
      }
    } else {
      k <- if (r == "4-1s") 4 else 10
      cut <- if (r == "4-1s") 1 else 0
      for (l in seq_len(n_levels)) {
        run <- which(roll_all(z[, l] > cut, k) | roll_all(z[, l] < -cut, k))
        add(run, r, rep(l, length(run)))
      }
    }
  }
  out <- if (length(rows) > 0) do.call(rbind, rows) else tibble::tibble(run = integer(), rule = character(), level = integer())
  out <- out[order(out$run, match(out$rule, rule_names), out$level), ]
  out$run <- as.integer(out$run)
  out$level <- as.integer(out$level)
  out
}
