# Run an expression with a given random seed and restore the previous state of the random number generator afterwards.
with_seed <- function(seed, code) {
  if (is.null(seed)) return(code)
  check_whole_seed <- is.numeric(seed) && length(seed) == 1 && is.finite(seed)
  if (!check_whole_seed) stop("`seed` must be NULL or a single number.", call. = FALSE)
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had) assign(".Random.seed", old, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
  }, add = TRUE)
  set.seed(seed)
  code
}

#' A randomized plan: completely randomized or randomized complete block
#'
#' Assigns treatments to runs **at random**. With `blocks = NULL` the design is **completely randomized**: all the runs are in one random order. With
#' `blocks` the design is a **randomized complete block design**: every block contains every treatment `replicates` times, and the order of the treatments
#' is randomized **within each block** (a block can be a boy with two feet, a day, a technician, a batch). Randomization protects against unknown
#' biases such as drift in time or a difference between left and right.
#'
#' @param treatments Treatment labels (at least two).
#' @param blocks `NULL` (completely randomized), a number of blocks, or a vector of block labels.
#' @param replicates Number of times each treatment appears in each block (or in the whole experiment if `blocks` is `NULL`). Default 1.
#' @param seed Optional seed, so that the plan can be reproduced. The state of the random number generator is restored afterwards.
#' @return A tibble with columns `run` (the order of the runs, in blocks if there are blocks), `block` (if any), and `treatment`.
#' @seealso [factorial_design()]
#' @examples
#' doe_plan(c("A", "B"), blocks = 5, seed = 1)       # the left or right foot of 5 boys, by coin flip
#' doe_plan(c("A", "B", "C"), replicates = 4, seed = 1)
#' @export
doe_plan <- function(treatments, blocks = NULL, replicates = 1, seed = NULL) {
  treatments <- as.character(treatments)
  if (length(unique(treatments)) < 2 || anyDuplicated(treatments)) stop("`treatments` must be at least two different labels.", call. = FALSE)
  check_count(replicates, "replicates", minimum = 1)
  with_seed(seed, {
    if (is.null(blocks)) {
      trt <- sample(rep(treatments, each = replicates))
      tibble::tibble(run = seq_along(trt), treatment = trt)
    } else {
      b <- if (length(blocks) == 1 && is.numeric(blocks)) { check_count(blocks, "blocks", minimum = 1); seq_len(blocks) } else blocks
      d <- do.call(rbind, lapply(b, function(bl) tibble::tibble(block = bl, treatment = sample(rep(treatments, each = replicates)))))
      d$run <- seq_len(nrow(d))
      d[c("run", "block", "treatment")]
    }
  })
}

# Pooled and paired t-test p-values for the rows of matrices (fast versions for the simulations)
t_pool_p <- function(x, y) {
  nx <- ncol(x); ny <- ncol(y)
  vx <- (rowSums(x^2) - nx * rowMeans(x)^2) / (nx - 1); vy <- (rowSums(y^2) - ny * rowMeans(y)^2) / (ny - 1)
  sp2 <- ((nx - 1) * vx + (ny - 1) * vy) / (nx + ny - 2)
  tt <- (rowMeans(x) - rowMeans(y)) / sqrt(sp2 * (1 / nx + 1 / ny))
  2 * stats::pt(-abs(tt), nx + ny - 2)
}
t_paired_p <- function(d) {
  n <- ncol(d); v <- (rowSums(d^2) - n * rowMeans(d)^2) / (n - 1)
  tt <- rowMeans(d) / sqrt(v / n)
  2 * stats::pt(-abs(tt), n - 1)
}

#' Power of the paired and the two-sample t-test as the units differ more
#'
#' Simulates the shoe experiment of Box, Hunter, and Hunter: each of `n` units (boys) tests **both** treatments (one on each foot). Unit `i` has its
#' own level \eqn{b_i \sim N(0, \sigma_b^2)} that is shared by its two results, the treatments differ by `delta`, and each result has an independent error
#' \eqn{N(0, \sigma_e^2)}. The same data are analyzed two ways: the **paired** t-test (the blocked analysis: the differences within units) and the
#' **two-sample** t-test (which ignores that each unit gave both results). Their power is estimated for each value of `sigma_b`.
#'
#' The paired test removes the unit-to-unit variation \eqn{\sigma_b}; the two-sample test does not, so its power falls as `sigma_b` grows while the power of
#' the paired test does not.
#'
#' @param n Number of units. Default 10.
#' @param delta Difference between the two treatments.
#' @param sigma_b Standard deviation(s) of the unit-to-unit variation (a vector gives one row per value).
#' @param sigma_e Standard deviation of the error of one result.
#' @param reps Number of simulated experiments for each `sigma_b`. Default 2000.
#' @param alpha Significance level. Default 0.05.
#' @return A tibble with `sigma_b`, `test` (`"paired"` or `"two-sample"`), and `power`.
#' @seealso [sim_randomization()], [sim_pseudorep()]
#' @examples
#' set.seed(1)
#' sim_blocking(n = 10, delta = 0.4, sigma_b = c(0.1, 1, 2), sigma_e = 0.3, reps = 500)
#' @export
sim_blocking <- function(n = 10, delta, sigma_b, sigma_e, reps = 2000, alpha = 0.05) {
  check_count(n, "n", minimum = 3); check_number(delta, "delta"); check_count(reps, "reps", minimum = 1); check_probability(alpha, "alpha")
  if (!(is.numeric(sigma_b) && all(is.finite(sigma_b)) && all(sigma_b >= 0))) stop("`sigma_b` must be non-negative numbers.", call. = FALSE)
  check_number(sigma_e, "sigma_e", positive = TRUE)
  out <- lapply(sigma_b, function(sb) {
    b <- matrix(stats::rnorm(reps * n, 0, sb), reps, n)
    a <- b + matrix(stats::rnorm(reps * n, 0, sigma_e), reps, n)
    bb <- b + delta + matrix(stats::rnorm(reps * n, 0, sigma_e), reps, n)
    tibble::tibble(sigma_b = sb, test = c("paired", "two-sample"), power = c(mean(t_paired_p(bb - a) < alpha), mean(t_pool_p(bb, a) < alpha)))
  })
  do.call(rbind, out)
}

#' The effect of randomization: a left-foot wear effect
#'
#' Simulates the shoe experiment with an extra **systematic** effect: the left foot wears the sole faster than the right by `foot_effect`. Each of `n` boys
#' tests material A on one foot and B on the other. With `randomize = TRUE` the foot of A is chosen by a coin flip for each boy; with `randomize = FALSE`
#' A is always on the left foot. The true difference of the materials is `delta` (B minus A). The paired t-test is used.
#'
#' Without randomization the foot effect is **confounded** with the material: the estimate is biased by `-foot_effect` and, when the true difference is 0,
#' the test finds a "significant" difference far more often than `alpha`. With randomization the foot effect is balanced and cancels on average.
#'
#' @param n Number of boys. Default 10.
#' @param delta True difference in wear, B minus A. Default 0.
#' @param foot_effect Extra wear of the left foot.
#' @param sigma_b,sigma_e Standard deviations of the boy-to-boy variation and of the error.
#' @param randomize Whether the foot of A is chosen at random. Default `TRUE`.
#' @param reps Number of simulated experiments. Default 2000.
#' @param alpha Significance level. Default 0.05.
#' @return A one-row tibble with `randomize`, `mean_estimate` (the average estimate of `delta`), `bias` (`mean_estimate - delta`), and `reject` (the share of
#'   experiments with a significant paired t-test).
#' @seealso [sim_blocking()], [doe_plan()]
#' @examples
#' set.seed(1)
#' sim_randomization(foot_effect = 0.5, sigma_b = 2, sigma_e = 0.3, randomize = FALSE, reps = 500)
#' sim_randomization(foot_effect = 0.5, sigma_b = 2, sigma_e = 0.3, randomize = TRUE, reps = 500)
#' @export
sim_randomization <- function(n = 10, delta = 0, foot_effect, sigma_b, sigma_e, randomize = TRUE, reps = 2000, alpha = 0.05) {
  check_count(n, "n", minimum = 3); check_number(delta, "delta"); check_number(foot_effect, "foot_effect"); check_number(sigma_b, "sigma_b")
  check_number(sigma_e, "sigma_e", positive = TRUE); check_count(reps, "reps", minimum = 1); check_probability(alpha, "alpha")
  if (!is.logical(randomize) || length(randomize) != 1 || is.na(randomize)) stop("`randomize` must be TRUE or FALSE.", call. = FALSE)
  b <- matrix(stats::rnorm(reps * n, 0, sigma_b), reps, n)
  a_left <- if (randomize) matrix(stats::rbinom(reps * n, 1, 0.5), reps, n) else matrix(1, reps, n)
  wear_a <- b + foot_effect * a_left + matrix(stats::rnorm(reps * n, 0, sigma_e), reps, n)
  wear_b <- b + delta + foot_effect * (1 - a_left) + matrix(stats::rnorm(reps * n, 0, sigma_e), reps, n)
  d <- wear_b - wear_a
  est <- mean(rowMeans(d))
  tibble::tibble(randomize = randomize, mean_estimate = est, bias = est - delta, reject = mean(t_paired_p(d) < alpha))
}

#' Pseudoreplication: duplicate measurements are not replicates
#'
#' Two treatments are compared with `n_units` units each (bags, patients, batches) and each unit is measured `m` times. The true difference is `delta`
#' (default 0). The simulation compares two analyses:
#' * **pseudoreplicated**: a two-sample t-test on all `n_units * m` measurements as if they were independent replicates;
#' * **correct**: a two-sample t-test on the `n_units` unit means.
#'
#' Repeated measurements of one unit estimate the **measurement** error, not the variation between units, so the pseudoreplicated test has far too many
#' degrees of freedom and, when `delta = 0`, a false alarm probability much larger than `alpha`.
#'
#' @param n_units Number of units per treatment. Default 5.
#' @param m Number of measurements of each unit. Default 3.
#' @param sigma_unit Standard deviation between units.
#' @param sigma_meas Standard deviation of the measurement error.
#' @param delta True difference between the treatments. Default 0.
#' @param reps Number of simulated experiments. Default 5000.
#' @param alpha Significance level. Default 0.05.
#' @return A tibble with `analysis` and `reject` (the share of experiments with a significant difference: the false alarm probability if `delta = 0`, the
#'   power otherwise).
#' @seealso [sim_blocking()]
#' @examples
#' set.seed(1)
#' sim_pseudorep(n_units = 5, m = 3, sigma_unit = 1, sigma_meas = 0.3, reps = 1000)
#' @export
sim_pseudorep <- function(n_units = 5, m = 3, sigma_unit, sigma_meas, delta = 0, reps = 5000, alpha = 0.05) {
  check_count(n_units, "n_units", minimum = 2); check_count(m, "m", minimum = 1); check_number(sigma_unit, "sigma_unit", positive = TRUE)
  check_number(sigma_meas, "sigma_meas", positive = TRUE); check_number(delta, "delta"); check_count(reps, "reps", minimum = 1); check_probability(alpha, "alpha")
  unit_a <- matrix(stats::rnorm(reps * n_units, 0, sigma_unit), reps, n_units)
  unit_b <- matrix(stats::rnorm(reps * n_units, delta, sigma_unit), reps, n_units)
  meas <- function(u) u[, rep(seq_len(n_units), each = m), drop = FALSE] + matrix(stats::rnorm(reps * n_units * m, 0, sigma_meas), reps, n_units * m)
  ma <- meas(unit_a); mb <- meas(unit_b)
  mean_by_unit <- function(x) { g <- rep(seq_len(n_units), each = m); sapply(seq_len(n_units), function(j) rowMeans(x[, g == j, drop = FALSE])) }
  ua <- mean_by_unit(ma); ub <- mean_by_unit(mb)
  if (n_units == 1) { ua <- matrix(ua, ncol = 1); ub <- matrix(ub, ncol = 1) }
  tibble::tibble(analysis = c("pseudoreplicated (all measurements)", "correct (unit means)"),
                 reject = c(mean(t_pool_p(ma, mb) < alpha), mean(t_pool_p(ua, ub) < alpha)))
}

#' A two-level factorial design
#'
#' Builds the run sheet of a **two-level factorial** experiment. Each factor has a low level (coded -1) and a high level (coded +1); the design has every
#' combination of the levels (the **standard order**: the first factor changes fastest), repeated `replicates` times, and the runs are put in a **random
#' order** (`run_order`) unless `randomize = FALSE`. With `generators` the design is a **fraction**: the generated factors are not crossed with the others
#' but defined as products of them (for example `list(C = c("A", "B"))` sets \eqn{C = AB} and gives the half fraction \eqn{2^{3-1}}).
#'
#' @param factors A character vector of factor names, or a named list whose elements are the two level labels of each factor, for example
#'   `list(wing = c("short", "long"), clip = c("no", "yes"))`.
#' @param replicates Number of replicates of every run. Default 1.
#' @param randomize Whether to randomize the run order. Default `TRUE`.
#' @param seed Optional seed for the randomization (the state of the random number generator is restored afterwards).
#' @param generators For a fraction: a named list; the name is a factor that is defined by the product of the factors in the element.
#' @return A tibble of class `factorial_design` with `run_order`, `std_order`, `replicate`, the factors coded -1 and +1, and (if labels were given) a
#'   column `<factor>_level` with the labels. Attributes `factors`, `levels`, and `generators` describe the design.
#' @seealso [factorial_effects()], [alias_structure()], [doe_plan()]
#' @examples
#' factorial_design(c("wing", "body", "clip"), replicates = 2, seed = 1)
#' factorial_design(c("A", "B", "C"), generators = list(C = c("A", "B")), randomize = FALSE)
#' @export
factorial_design <- function(factors, replicates = 1, randomize = TRUE, seed = NULL, generators = NULL) {
  labels <- NULL
  if (is.list(factors)) {
    labels <- factors
    if (is.null(names(labels)) || any(lengths(labels) != 2)) stop("`factors` must be names, or a named list of two level labels per factor.", call. = FALSE)
    factors <- names(labels)
  }
  if (!is.character(factors) || length(factors) < 2 || anyDuplicated(factors) || anyNA(factors)) stop("`factors` must be at least two different names.", call. = FALSE)
  check_count(replicates, "replicates", minimum = 1)
  if (!is.logical(randomize) || length(randomize) != 1 || is.na(randomize)) stop("`randomize` must be TRUE or FALSE.", call. = FALSE)
  gen_names <- names(generators)
  if (length(generators)) {
    if (is.null(gen_names) || !all(gen_names %in% factors)) stop("The names of `generators` must be factors.", call. = FALSE)
    for (g in generators) if (!all(g %in% setdiff(factors, gen_names)) || length(g) < 2) stop("Each generator must be a product of at least two base factors.", call. = FALSE)
  }
  base <- setdiff(factors, gen_names)
  kb <- length(base)
  grid <- as.matrix(expand.grid(rep(list(c(-1, 1)), kb)))
  colnames(grid) <- base
  des <- as.data.frame(grid)
  for (g in gen_names) des[[g]] <- Reduce(`*`, des[generators[[g]]])
  des <- des[factors]
  des$std_order <- seq_len(nrow(des))
  des <- des[rep(seq_len(nrow(des)), times = replicates), , drop = FALSE]
  des$replicate <- rep(seq_len(replicates), each = nrow(des) / replicates)
  des <- with_seed(seed, if (randomize) des[sample(nrow(des)), , drop = FALSE] else des)
  des$run_order <- seq_len(nrow(des))
  out <- des[c("run_order", "std_order", "replicate", factors)]
  if (!is.null(labels)) for (f in factors) out[[paste0(f, "_level")]] <- ifelse(out[[f]] < 0, labels[[f]][1], labels[[f]][2])
  out <- tibble::as_tibble(out)
  attr(out, "factors") <- factors; attr(out, "levels") <- labels; attr(out, "generators") <- generators
  class(out) <- c("factorial_design", class(out))
  out
}

#' @export
print.factorial_design <- function(x, ...) {
  f <- attr(x, "factors"); g <- attr(x, "generators")
  cat("<factorial_design> ", length(f), " factors (", paste(f, collapse = ", "), "), ", nrow(x), " runs",
      if (length(g)) paste0(", fraction with ", paste(names(g), "=", vapply(g, paste, "", collapse = ""), collapse = ", ")), "\n", sep = "")
  NextMethod()
}

# Effects (sets of factor names) of order up to `max_order`, as character vectors
effect_sets <- function(factors, max_order) {
  out <- list()
  for (o in seq_len(min(max_order, length(factors)))) out <- c(out, utils::combn(factors, o, simplify = FALSE))
  out
}
effect_label <- function(s) paste(s, collapse = ":")

#' The alias structure of a fractional factorial design
#'
#' In a fraction the effects cannot all be estimated: some are **aliased** (confounded), so that their estimates are the sum of the effects. The
#' **defining relation** is the set of products of the generators that equal the identity; multiplying an effect by each word of the defining relation
#' gives the effects it is aliased with. For \eqn{C = AB} the defining relation is \eqn{I = ABC}, so \eqn{A = BC}, \eqn{B = AC}, and \eqn{C = AB}.
#'
#' @param factors Character vector of the factor names.
#' @param generators A named list as in [factorial_design()].
#' @return A tibble with one row per group of aliased effects: `aliases` (the effects of the group, joined by `=`) and `order` (the smallest order of
#'   an effect in the group). The attribute `resolution` is the length of the shortest word of the defining relation, and `defining` lists the words.
#' @seealso [factorial_design()]
#' @examples
#' alias_structure(c("A", "B", "C"), list(C = c("A", "B")))
#' alias_structure(c("A", "B", "C", "D"), list(D = c("A", "B", "C")))
#' @export
alias_structure <- function(factors, generators) {
  if (!is.character(factors) || length(factors) < 3) stop("`factors` must be at least three names.", call. = FALSE)
  if (!is.list(generators) || length(generators) == 0 || is.null(names(generators)) || !all(names(generators) %in% factors)) stop("`generators` must be a named list of factors.", call. = FALSE)
  words <- lapply(names(generators), function(g) sort(c(g, generators[[g]])))
  sym <- function(a, b) sort(union(setdiff(a, b), setdiff(b, a)))
  # all products of subsets of the generator words
  defining <- list()
  for (i in seq_along(words)) for (s in utils::combn(seq_along(words), i, simplify = FALSE)) defining <- c(defining, list(Reduce(sym, words[s])))
  defining <- defining[!duplicated(lapply(defining, paste, collapse = "|"))]
  lens <- vapply(defining, length, 1L)
  effects <- effect_sets(factors, length(factors) - 1L)
  key <- function(s) paste(sort(s), collapse = "|")
  groups <- list(); seen <- character()
  for (e in effects) {
    if (key(e) %in% seen) next
    grp <- c(list(sort(e)), lapply(defining, function(w) sym(e, w)))
    grp <- grp[lengths(grp) > 0]
    grp <- grp[!duplicated(vapply(grp, key, ""))]
    seen <- c(seen, vapply(grp, key, ""))
    grp <- grp[order(lengths(grp), vapply(grp, key, ""))]
    groups <- c(groups, list(grp))
  }
  out <- tibble::tibble(aliases = vapply(groups, function(g) paste(vapply(g, effect_label, ""), collapse = " = "), ""), order = vapply(groups, function(g) min(lengths(g)), 1L))
  out <- out[order(out$order), ]
  attr(out, "resolution") <- min(lens)
  attr(out, "defining") <- vapply(defining, function(w) paste(w, collapse = ""), "")
  out
}

# Code a column as -1 / +1
code_pm1 <- function(v, name) {
  if (is.factor(v) || is.character(v)) {
    lv <- if (is.factor(v)) levels(v) else sort(unique(v))
    if (length(unique(v)) != 2) stop("Factor `", name, "` must have exactly two levels.", call. = FALSE)
    return(ifelse(v == lv[1], -1, 1))
  }
  u <- sort(unique(v))
  if (length(u) != 2) stop("Factor `", name, "` must have exactly two levels.", call. = FALSE)
  ifelse(v == u[1], -1, 1)
}

#' Effects of a two-level factorial experiment
#'
#' Estimates the **main effects and interactions** of a two-level factorial experiment by contrasts. With the factors coded -1 and +1, the **effect** of a
#' factor is the average response at its high level minus the average at its low level, which is twice the coefficient of the regression of the response on
#' the coded factor. An **interaction** effect is the effect of the product of the coded factors: for two factors, half the difference between the effect
#' of A at the high level of B and its effect at the low level of B. For a **replicated** experiment the error standard deviation (pure error) estimates
#' the standard error of every effect, \eqn{2\sigma/\sqrt{N}}, and gives t-tests and confidence intervals. For an **unreplicated** saturated experiment
#' there is no error estimate and the effects are judged by the half-normal plot of [autoplot.factorial_effects()].
#'
#' @param data A data frame with the response and the factors.
#' @param response The response column (unquoted).
#' @param factors Character vector of the factor columns. Two-level numeric, character, or factor columns are coded -1 (the lower or first level) and +1.
#' @param order Highest order of the interactions in the model. Default: all.
#' @return An object of class `factorial_effects`: a list with `effects` (a tibble with `term`, `effect`, `se`, `t`, `p`, `lcl`, `ucl`), `coefficients`
#'   (the regression coefficients of the coded model, half the effects), `intercept` (the grand mean), `sigma` (the error standard deviation; `NA` without
#'   replication), `df` (its degrees of freedom), `n`, `factors`, `terms`, and `aliased` (terms that cannot be estimated, in a fraction). Use
#'   [predict()] for predictions and [autoplot()] for the plots.
#' @seealso [factorial_design()], [interaction_plot()], [doe_power()]
#' @examples
#' fe <- factorial_effects(helicopter, flight_time, c("wing", "body", "clip"))
#' fe
#' predict(fe, data.frame(wing = 1, body = 1, clip = 1))
#' @export
factorial_effects <- function(data, response, factors, order = NULL) {
  data <- as.data.frame(data)
  y <- rlang::eval_tidy(rlang::enquo(response), data)
  if (!is.numeric(y) || anyNA(y)) stop("`response` must be a numeric column with no missing values.", call. = FALSE)
  if (!is.character(factors) || length(factors) < 2 || !all(factors %in% names(data))) stop("`factors` must name at least two columns of `data`.", call. = FALSE)
  X <- vapply(factors, function(f) code_pm1(data[[f]], f), numeric(nrow(data)))
  if (is.null(dim(X))) X <- matrix(X, ncol = length(factors), dimnames = list(NULL, factors))
  max_order <- if (is.null(order)) length(factors) else { check_count(order, "order", minimum = 1); order }
  sets <- effect_sets(factors, max_order)
  terms <- vapply(sets, effect_label, "")
  M <- vapply(sets, function(s) Reduce(`*`, lapply(s, function(f) X[, f])), numeric(nrow(X)))
  colnames(M) <- terms
  fit <- stats::lm.fit(cbind(`(Intercept)` = 1, M), y)
  ok <- !is.na(fit$coefficients)
  est <- terms[ok[-1]]; aliased <- terms[!ok[-1]]
  Mk <- cbind(`(Intercept)` = 1, M[, est, drop = FALSE])
  fit <- stats::lm.fit(Mk, y)
  n <- length(y); p <- ncol(Mk); df <- n - p
  coefs <- fit$coefficients
  sigma <- if (df > 0) sqrt(sum(fit$residuals^2) / df) else NA_real_
  se_coef <- if (df > 0) sigma * sqrt(diag(solve(crossprod(Mk)))) else rep(NA_real_, p)
  eff <- 2 * coefs[-1]; se <- 2 * se_coef[-1]
  tt <- eff / se; tcrit <- if (df > 0) stats::qt(0.975, df) else NA_real_
  tab <- tibble::tibble(term = est, effect = unname(eff), se = unname(se), t = unname(tt), p = if (df > 0) unname(2 * stats::pt(-abs(tt), df)) else NA_real_,
                        lcl = unname(eff - tcrit * se), ucl = unname(eff + tcrit * se))
  structure(list(effects = tab, coefficients = stats::setNames(unname(coefs[-1]), est), intercept = unname(coefs[1]), sigma = sigma, df = df, n = n,
                 factors = factors, terms = est, aliased = aliased), class = "factorial_effects")
}

#' @export
print.factorial_effects <- function(x, ...) {
  cat("<factorial_effects> ", x$n, " runs, ", length(x$factors), " factors; grand mean ", format(signif(x$intercept, 5)), "\n", sep = "")
  if (is.na(x$sigma)) cat("  no error estimate (no replication, saturated model): use the half-normal plot\n") else
    cat("  error standard deviation ", format(signif(x$sigma, 3)), " on ", x$df, " degrees of freedom; standard error of an effect ", format(signif(2 * x$sigma / sqrt(x$n), 3)), "\n", sep = "")
  tab <- as.data.frame(x$effects); num <- vapply(tab, is.numeric, TRUE); tab[num] <- lapply(tab[num], function(v) signif(v, 3))
  print(tab, row.names = FALSE)
  if (length(x$aliased)) cat("  not estimable (aliased): ", paste(x$aliased, collapse = ", "), "\n", sep = "")
  invisible(x)
}

#' @export
as.data.frame.factorial_effects <- function(x, ...) as.data.frame(x$effects)

#' Predictions of a factorial model
#'
#' @param object A [factorial_effects()] object.
#' @param newdata A data frame with the factors coded -1 and +1 (the names of the factors).
#' @param ... Not used.
#' @return A numeric vector of predicted responses: the grand mean plus half of each effect times the product of the coded factors.
#' @importFrom stats predict
#' @export
predict.factorial_effects <- function(object, newdata, ...) {
  newdata <- as.data.frame(newdata)
  if (!all(object$factors %in% names(newdata))) stop("`newdata` must have the factors ", paste(object$factors, collapse = ", "), ".", call. = FALSE)
  if (!all(vapply(object$factors, function(f) all(newdata[[f]] %in% c(-1, 1)), TRUE))) stop("The factors in `newdata` must be coded -1 and +1.", call. = FALSE)
  pred <- rep(object$intercept, nrow(newdata))
  for (tm in object$terms) pred <- pred + object$coefficients[[tm]] * Reduce(`*`, lapply(strsplit(tm, ":", fixed = TRUE)[[1]], function(f) newdata[[f]]))
  pred
}

#' Draw the effects of a factorial experiment with ggplot2
#'
#' `which = "effects"` draws every effect with its 95% confidence interval (for a replicated experiment), ordered by size, around the line at 0.
#' `which = "half_normal"` draws the absolute effects against the half-normal quantiles: effects that are only noise fall along a line through the
#' origin and the **active** effects are far above it.
#'
#' @param object A [factorial_effects()] object.
#' @param which `"effects"` (default) or `"half_normal"`.
#' @param title,subtitle Plot title and subtitle.
#' @param ... Not used.
#' @return A ggplot object.
#' @export
autoplot.factorial_effects <- function(object, which = c("effects", "half_normal"), title = NULL, subtitle = NULL, ...) {
  which <- match.arg(which)
  d <- object$effects
  if (which == "effects") {
    if (is.na(object$sigma)) stop("The effects plot needs a replicated experiment; use `which = \"half_normal\"`.", call. = FALSE)
    d$term <- factor(d$term, levels = d$term[order(abs(d$effect))])
    g <- ggplot2::ggplot(d, ggplot2::aes(x = .data$effect, y = .data$term)) + ggplot2::geom_vline(xintercept = 0) +
      ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$lcl, xmax = .data$ucl), width = 0.25, orientation = "y", color = "#4c78a8") +
      ggplot2::geom_point(ggplot2::aes(color = .data$p < 0.05), size = 2.2) + ggplot2::scale_color_manual(values = c(`TRUE` = "#e45756", `FALSE` = "black"), guide = "none") +
      ggplot2::labs(x = "Effect (high level minus low level) with 95% confidence interval", y = NULL, title = if (is.null(title)) "Effects" else title, subtitle = subtitle)
  } else {
    a <- sort(abs(d$effect)); nn <- length(a)
    q <- stats::qnorm(0.5 + 0.5 * (seq_len(nn) - 0.5) / nn)
    lab <- d$term[order(abs(d$effect))]
    dd <- tibble::tibble(q = q, a = a, lab = lab)
    slope <- sum(dd$q[seq_len(max(1, floor(nn * 0.6)))] * dd$a[seq_len(max(1, floor(nn * 0.6)))]) / sum(dd$q[seq_len(max(1, floor(nn * 0.6)))]^2)
    g <- ggplot2::ggplot(dd, ggplot2::aes(x = .data$q, y = .data$a)) + ggplot2::geom_abline(slope = slope, intercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::geom_point(size = 2) + ggplot2::geom_text(ggplot2::aes(label = .data$lab), hjust = -0.15, vjust = 0.4, size = 3.2) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.25))) +
      ggplot2::labs(x = "Half-normal quantile", y = "Absolute effect", title = if (is.null(title)) "Half-normal plot of the effects" else title, subtitle = subtitle)
  }
  g + ggplot2::theme_bw() + ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

#' Interaction plot
#'
#' Plots the **mean response** at each level of `x`, with one line for each level of `trace`. **Parallel lines** mean no interaction (the effect of `x`
#' is the same at every level of `trace`); lines that are not parallel, or that cross, mean an interaction. The individual observations are shown as
#' light points.
#'
#' @param data A data frame.
#' @param response The response column (unquoted).
#' @param x,trace Names of the two factor columns (quoted).
#' @param title,subtitle,xlab,ylab Plot title, subtitle, and axis labels.
#' @return A ggplot object.
#' @seealso [factorial_effects()]
#' @examples
#' interaction_plot(helicopter, flight_time, "wing", "clip")
#' @export
interaction_plot <- function(data, response, x, trace, title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL) {
  data <- as.data.frame(data)
  y <- rlang::eval_tidy(rlang::enquo(response), data)
  if (!is.character(x) || !is.character(trace) || length(x) != 1 || length(trace) != 1 || !all(c(x, trace) %in% names(data))) stop("`x` and `trace` must name columns of `data`.", call. = FALSE)
  d <- data.frame(x = factor(data[[x]]), trace = factor(data[[trace]]), y = y)
  m <- stats::aggregate(y ~ x + trace, d, mean)
  ggplot2::ggplot(m, ggplot2::aes(x = .data$x, y = .data$y, color = .data$trace, group = .data$trace)) +
    ggplot2::geom_point(data = d, ggplot2::aes(x = .data$x, y = .data$y, color = .data$trace), alpha = 0.35, size = 1.4, position = ggplot2::position_dodge(width = 0.12), inherit.aes = FALSE) +
    ggplot2::geom_line(linewidth = 1) + ggplot2::geom_point(size = 3) +
    ggplot2::labs(x = if (is.null(xlab)) x else xlab, y = if (is.null(ylab)) "Mean response" else ylab, color = trace, title = title, subtitle = subtitle) +
    ggplot2::theme_bw() + ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

#' Power of an effect test in a two-level factorial experiment
#'
#' With `N` runs in a two-level factorial in `k` factors (all the interactions in the model, so `N - 2^k` degrees of freedom for error), the standard error of
#' an effect is \eqn{2\sigma/\sqrt{N}}. The power of the two-sided t-test of an effect of size `delta` at level `alpha` follows from the noncentral t
#' distribution with noncentrality \eqn{\delta\sqrt{N}/(2\sigma)}. `doe_n()` tabulates the power for 2, 3, ... replicates.
#'
#' @param delta The effect to detect (high level minus low level).
#' @param sigma Standard deviation of the error of one run.
#' @param N Total number of runs (a multiple of \eqn{2^k}).
#' @param k Number of factors. Default 3.
#' @param alpha Significance level. Default 0.05.
#' @param power Target power for `doe_n()`. Default 0.8.
#' @param max_replicates Largest number of replicates in the table of `doe_n()`. Default 12.
#' @return `doe_power()`: the power. `doe_n()`: a tibble with `replicates`, `N`, `se` (standard error of an effect), `power`, and `ok` (power at least `power`).
#' @seealso [factorial_effects()]
#' @examples
#' doe_power(delta = 2, sigma = 2, N = 32)
#' doe_n(delta = 2, sigma = 2)
#' @export
doe_power <- function(delta, sigma, N, k = 3, alpha = 0.05) {
  check_number(delta, "delta"); check_number(sigma, "sigma", positive = TRUE); check_count(k, "k", minimum = 1); check_probability(alpha, "alpha")
  if (!(is.numeric(N) && all(is.finite(N)) && all(N > 2^k))) stop("`N` must be larger than 2^k, the number of parameters.", call. = FALSE)
  df <- N - 2^k
  ncp <- abs(delta) / (2 * sigma / sqrt(N)); tc <- stats::qt(1 - alpha / 2, df)
  1 - stats::pt(tc, df, ncp) + stats::pt(-tc, df, ncp)
}

#' @rdname doe_power
#' @export
doe_n <- function(delta, sigma, power = 0.8, k = 3, alpha = 0.05, max_replicates = 12) {
  check_probability(power, "power"); check_count(max_replicates, "max_replicates", minimum = 2)
  target <- power
  r <- 2:max_replicates; N <- r * 2^k
  pw <- doe_power(delta, sigma, N, k = k, alpha = alpha)
  tibble::tibble(replicates = r, N = N, se = 2 * sigma / sqrt(N), power = pw, ok = pw >= target)
}

#' One-factor-at-a-time against a factorial design: a simulation
#'
#' Two factors at two levels (-1, +1) have a true response (to be maximized) with main effects `b1` and `b2` and an interaction `b12`: the mean at
#' \eqn{(x_1, x_2)} is \eqn{(b_1x_1 + b_2x_2 + b_{12}x_1x_2)/2} (so `b1` is the effect of factor 1 when the interaction is 0). Each combination that is run
#' is run `r` times, with error standard deviation `sigma`. Two strategies pick the best combination:
#' * **one factor at a time**: start at (-1, -1); change factor 1 and keep the better level; then change factor 2 at that level and keep the better level
#'   (three combinations are run);
#' * **factorial**: run all four combinations and pick the one with the highest mean.
#'
#' With a large interaction the effect of a factor at the starting level of the other can be small or the wrong sign, so the one-factor-at-a-time path
#' never visits the best combination.
#'
#' @param b1,b2,b12 True effects of factor 1, factor 2, and their interaction.
#' @param sigma Standard deviation of the error of one run.
#' @param r Number of runs at each combination. Default 1.
#' @param reps Number of simulated experiments. Default 5000.
#' @return A tibble with `strategy`, `runs` (the number of runs), `p_best` (the probability of choosing the truly best combination), and `loss` (the
#'   average shortfall of the chosen combination from the best, in response units). The attribute `means` has the true means of the four combinations.
#' @seealso [factorial_effects()]
#' @examples
#' set.seed(1)
#' sim_ofat(b1 = 0.3, b2 = -0.1, b12 = 0.4, sigma = 0.15, r = 2, reps = 1000)
#' @export
sim_ofat <- function(b1, b2, b12, sigma, r = 1, reps = 5000) {
  check_number(b1, "b1"); check_number(b2, "b2"); check_number(b12, "b12"); check_number(sigma, "sigma", positive = TRUE)
  check_count(r, "r", minimum = 1); check_count(reps, "reps", minimum = 1)
  a <- c(-1, 1, -1, 1); cc <- c(-1, -1, 1, 1)              # combinations: (-,-), (+,-), (-,+), (+,+)
  mu <- (b1 * a + b2 * cc + b12 * a * cc) / 2
  best <- which.max(mu)
  if (sum(mu == max(mu)) > 1) warning("Two combinations are equally good.", call. = FALSE)
  obs <- matrix(mu, reps, 4, byrow = TRUE) + matrix(stats::rnorm(reps * 4, 0, sigma / sqrt(r)), reps, 4)
  # one factor at a time: start at (-,-)=1, try (+,-)=2, then change factor 2 at the chosen level of factor 1
  step1 <- ifelse(obs[, 2] > obs[, 1], 2, 1)
  alt <- ifelse(step1 == 2, 4, 3)
  idx <- cbind(seq_len(reps), step1); idx2 <- cbind(seq_len(reps), alt)
  ofat <- ifelse(obs[idx2] > obs[idx], alt, step1)
  fact <- max.col(obs, ties.method = "first")
  res <- tibble::tibble(strategy = c("one factor at a time", "factorial"), runs = c(3 * r, 4 * r),
                        p_best = c(mean(ofat == best), mean(fact == best)), loss = c(mean(max(mu) - mu[ofat]), mean(max(mu) - mu[fact])))
  attr(res, "means") <- stats::setNames(mu, c("(-,-)", "(+,-)", "(-,+)", "(+,+)"))
  res
}

#' Cause-and-effect (fishbone) diagram
#'
#' Draws a **cause-and-effect diagram** (Ishikawa or fishbone): the effect (the problem or the response) at the head, the main categories of causes as
#' bones, and the individual causes along each bone. It is a brainstorming tool to organize the factors that may influence a response before choosing
#' the ones to put in an experiment.
#'
#' @param effect The effect (problem or response), a short text.
#' @param causes A named list: each name is a category (for example `"Materials"`) and each element a character vector of causes.
#' @param title Plot title.
#' @return A ggplot object.
#' @examples
#' fishbone("Short flight time",
#'          list(Materials = c("paper weight", "clip"),
#'               Design = c("wing length", "body length", "wing width"),
#'               Method = c("drop height", "release"),
#'               Measurement = c("stopwatch", "observer")))
#' @export
fishbone <- function(effect, causes, title = NULL) {
  if (!is.character(effect) || length(effect) != 1) stop("`effect` must be a single text.", call. = FALSE)
  if (!is.list(causes) || length(causes) < 2 || is.null(names(causes))) stop("`causes` must be a named list with at least two categories.", call. = FALSE)
  nc <- length(causes); n_top <- ceiling(nc / 2)
  spine_end <- 10 + 3 * max(0, n_top - 3)
  xs <- seq(spine_end * 0.28, spine_end * 0.88, length.out = n_top)
  bones <- list(); lines <- list(); labs <- list(); cats <- list()
  for (i in seq_len(nc)) {
    top <- i %% 2 == 1; j <- ceiling(i / 2); xe <- xs[j]; sgn <- if (top) 1 else -1
    x0 <- xe - 2.4; y0 <- sgn * 4.2
    bones[[i]] <- data.frame(x = x0, y = y0, xend = xe, yend = 0)
    cats[[i]] <- data.frame(x = x0, y = y0 + sgn * 0.45, label = names(causes)[i])
    cs <- causes[[i]]; m <- length(cs)
    if (m > 0) for (q in seq_len(m)) {
      tq <- q / (m + 1)
      bx <- x0 + (xe - x0) * tq; by <- y0 * (1 - tq)
      lines[[length(lines) + 1]] <- data.frame(x = bx - 1.6, y = by, xend = bx, yend = by)
      labs[[length(labs) + 1]] <- data.frame(x = bx - 1.65, y = by, label = cs[q])
    }
  }
  bones <- do.call(rbind, bones); cats <- do.call(rbind, cats); lines <- do.call(rbind, lines); labs <- do.call(rbind, labs)
  ggplot2::ggplot() +
    ggplot2::geom_segment(ggplot2::aes(x = 0, y = 0, xend = spine_end, yend = 0), linewidth = 1.2, arrow = ggplot2::arrow(length = ggplot2::unit(0.15, "inches"), type = "closed")) +
    ggplot2::geom_segment(data = bones, ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend), linewidth = 0.9, color = "#4c78a8") +
    ggplot2::geom_segment(data = lines, ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend), linewidth = 0.4, color = "grey40") +
    ggplot2::geom_text(data = labs, ggplot2::aes(x = .data$x, y = .data$y, label = .data$label), hjust = 1, vjust = -0.35, size = 3.3) +
    ggplot2::geom_label(data = cats, ggplot2::aes(x = .data$x, y = .data$y, label = .data$label), fontface = "bold", fill = "#e8eef6", linewidth = 0.3) +
    ggplot2::annotate("label", x = spine_end + 0.2, y = 0, label = effect, hjust = 0, fontface = "bold", fill = "#fdebd5", size = 4) +
    ggplot2::coord_cartesian(xlim = c(-3.5, spine_end + 0.3 + 0.28 * nchar(effect)), ylim = c(-6, 6), clip = "off") +
    ggplot2::labs(title = title) + ggplot2::theme_void() + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
}
