#' Single sampling plan
#'
#' A single sampling plan for acceptance sampling by attributes: inspect `n` items drawn at
#' random from a lot, and **accept** the lot if the number of defective items found is at most
#' `c` (the acceptance number); otherwise **reject** it. If the lot size `N` is given, the
#' number of defectives in the sample is hypergeometric; if not, the binomial approximation
#' (a very large lot, or sampling with replacement) is used.
#'
#' @param n Sample size (whole number, at least 1).
#' @param c Acceptance number (whole number, at least 0 and less than `n`).
#' @param N Lot size (whole number at least `n`), or `NULL` (the default) for the binomial
#'   approximation.
#' @return An object of class `sampling_plan` with elements `n`, `c`, and `N`.
#' @examples
#' plan <- sampling_plan(n = 50, c = 1)
#' plan
#' sampling_plan(n = 50, c = 1, N = 200)
#' @seealso [oc_curve()], [plan_risks()], [find_plan()]
#' @export
sampling_plan <- function(n, c, N = NULL) {
  check_whole(n, "n", 1)
  check_whole(c, "c", 0)
  if (c >= n) stop("`c` must be smaller than `n`; otherwise every lot is accepted.", call. = FALSE)
  if (!is.null(N)) {
    check_whole(N, "N", n)
  }
  structure(list(n = n, c = c, N = N), class = "sampling_plan")
}

#' @export
print.sampling_plan <- function(x, ...) {
  cat("<sampling_plan> inspect n = ", x$n, " items per lot; accept the lot if at most c = ", x$c,
      " defective", if (x$c != 1) "s" else "", " are found\n", sep = "")
  if (is.null(x$N)) {
    cat("  number of defectives in the sample: binomial (large lot)\n")
  } else {
    cat("  lot size N = ", x$N, "; number of defectives in the sample: hypergeometric\n", sep = "")
  }
  invisible(x)
}

# Probability of acceptance of a plan at lot fractions defective `p`.
plan_pa <- function(plan, p) {
  if (!is.numeric(p) || anyNA(p) || any(p < 0) || any(p > 1)) {
    stop("`p` must be fractions defective between 0 and 1.", call. = FALSE)
  }
  if (is.null(plan$N)) {
    stats::pbinom(plan$c, plan$n, p)
  } else {
    d <- round(plan$N * p)        # number of defectives in the lot
    stats::phyper(plan$c, d, plan$N - d, plan$n)
  }
}

check_plan <- function(plan) {
  if (!inherits(plan, "sampling_plan")) stop("`plan` must be created with sampling_plan().", call. = FALSE)
  invisible(plan)
}

#' Operating characteristic (OC) curve of a sampling plan
#'
#' The OC curve gives the probability of acceptance \eqn{P_a(p)}{Pa(p)} of a lot as a function
#' of the lot's true fraction defective \eqn{p}{p}. For a plan that accepts when the sample has at
#' most `c` defectives,
#' \deqn{P_a(p) = \sum_{x=0}^{c} \binom{n}{x} p^x (1-p)^{n-x}}{Pa(p) = sum_{x = 0}^{c} choose(n, x) p^x (1 - p)^(n - x)}
#' (binomial), or the hypergeometric analogue when the lot size is finite. With a lot size `N`,
#' the number of defectives in the lot is `round(N * p)`.
#'
#' @param plan A [sampling_plan()], or a named list of plans to compare.
#' @param p Fractions defective at which to evaluate the curve. Default `seq(0, 0.2, by = 0.002)`.
#' @return A tibble of class `oc_curve` with columns `plan` (a label), `p`, and `pa`.
#' @examples
#' oc <- oc_curve(sampling_plan(n = 50, c = 1))
#' head(oc)
#' oc_curve(list("n = 50, c = 1" = sampling_plan(50, 1), "n = 100, c = 2" = sampling_plan(100, 2)),
#'          p = c(0.01, 0.05))
#' @seealso [autoplot.oc_curve()], [plan_risks()]
#' @export
oc_curve <- function(plan, p = seq(0, 0.2, by = 0.002)) {
  plans <- if (inherits(plan, "sampling_plan")) list(plan) else plan
  if (!is.list(plans) || length(plans) == 0) stop("`plan` must be a sampling_plan() or a list of them.", call. = FALSE)
  labels <- names(plans)
  if (is.null(labels)) labels <- rep("", length(plans))
  out <- lapply(seq_along(plans), function(i) {
    check_plan(plans[[i]])
    label <- if (labels[i] != "") labels[i] else
      paste0("n = ", plans[[i]]$n, ", c = ", plans[[i]]$c)
    tibble::tibble(plan = label, p = p, pa = plan_pa(plans[[i]], p))
  })
  res <- do.call(rbind, out)
  res$plan <- factor(res$plan, levels = unique(res$plan))
  structure(res, class = c("oc_curve", class(res)))
}

#' Draw an OC curve with ggplot2
#'
#' @param object An [oc_curve()] object.
#' @param aql,ltpd Optional acceptable quality level and lot tolerance percent defective (as
#'   fractions). If given, the curve is marked at these points and the producer's and consumer's
#'   risks are shown (for the first plan only).
#' @param title,subtitle,xlab,ylab Plot title, subtitle, and axis labels.
#' @param ... Not used.
#' @return A ggplot object.
#' @examples
#' plan <- sampling_plan(50, 1)
#' autoplot(oc_curve(plan), aql = 0.01, ltpd = 0.05)
#' @export
autoplot.oc_curve <- function(object, aql = NULL, ltpd = NULL, title = "OC curve",
                              subtitle = NULL, xlab = "Lot fraction defective, p",
                              ylab = "Probability of accepting the lot", ...) {
  pct <- function(x) paste0(format(100 * x, trim = TRUE), "%")
  df <- as.data.frame(object)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$p, y = .data$pa, color = .data$plan)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_x_continuous(labels = pct) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
    ggplot2::labs(x = xlab, y = ylab, title = title, subtitle = subtitle, color = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = if (nlevels(df$plan) > 1) "bottom" else "none")
  first <- df[df$plan == levels(df$plan)[1], ]
  mark <- function(p, pa, label) {
    list(ggplot2::annotate("segment", x = p, xend = p, y = 0, yend = pa, linetype = "dotted"),
         ggplot2::annotate("segment", x = 0, xend = p, y = pa, yend = pa, linetype = "dotted"),
         ggplot2::annotate("point", x = p, y = pa, size = 2.6),
         ggplot2::annotate("text", x = p, y = pa, label = label, hjust = -0.1, vjust = -0.6, size = 3.3))
  }
  if (!is.null(aql)) {
    pa <- stats::approx(first$p, first$pa, xout = aql)$y
    p <- p + mark(aql, pa, sprintf("AQL: Pa = %.3f (producer's risk %.3f)", pa, 1 - pa))
  }
  if (!is.null(ltpd)) {
    pa <- stats::approx(first$p, first$pa, xout = ltpd)$y
    p <- p + mark(ltpd, pa, sprintf("LTPD: Pa = %.3f (consumer's risk)", pa))
  }
  p
}

#' Producer's and consumer's risk of a sampling plan
#'
#' The **producer's risk** \eqn{\alpha}{alpha} is the probability that a lot of acceptable quality
#' (fraction defective equal to the acceptable quality level, AQL) is rejected,
#' \eqn{\alpha = 1 - P_a(\text{AQL})}{alpha = 1 - Pa(AQL)}. The **consumer's risk**
#' \eqn{\beta}{beta} is the probability that a lot of unacceptable quality (the lot tolerance
#' percent defective, LTPD) is accepted, \eqn{\beta = P_a(\text{LTPD})}{beta = Pa(LTPD)}.
#'
#' @param plan A [sampling_plan()].
#' @param aql Acceptable quality level, as a fraction defective.
#' @param ltpd Lot tolerance percent defective, as a fraction (greater than `aql`).
#' @return A one-row tibble with `aql`, `ltpd`, `pa_aql`, `pa_ltpd`, `producers_risk`, and
#'   `consumers_risk`.
#' @examples
#' plan_risks(sampling_plan(50, 1), aql = 0.01, ltpd = 0.05)
#' @seealso [oc_curve()], [find_plan()]
#' @export
plan_risks <- function(plan, aql, ltpd) {
  check_plan(plan)
  if (!is.numeric(aql) || !is.numeric(ltpd) || length(aql) != 1 || length(ltpd) != 1 ||
      aql <= 0 || ltpd >= 1 || aql >= ltpd) {
    stop("`aql` and `ltpd` must be single fractions with 0 < aql < ltpd < 1.", call. = FALSE)
  }
  pa_aql <- plan_pa(plan, aql)
  pa_ltpd <- plan_pa(plan, ltpd)
  tibble::tibble(aql = aql, ltpd = ltpd, pa_aql = pa_aql, pa_ltpd = pa_ltpd,
                 producers_risk = 1 - pa_aql, consumers_risk = pa_ltpd)
}

check_rectifying <- function(plan) {
  check_plan(plan)
  if (is.null(plan$N)) stop("Rectifying inspection needs the lot size: use sampling_plan(n, c, N = ...).", call. = FALSE)
  invisible(plan)
}

#' Average outgoing quality and average total inspection (rectifying inspection)
#'
#' With **rectifying inspection**, every rejected lot is inspected 100 percent and every
#' defective item found is replaced by a good one; accepted lots are shipped as they are (after
#' the defectives found in the sample are replaced). For a lot of size `N` with fraction defective `p`:
#' \deqn{\text{AOQ}(p) = P_a(p)\,p\,\frac{N-n}{N}, \qquad \text{ATI}(p) = n + \bigl(1 - P_a(p)\bigr)(N - n).}{AOQ(p) = Pa(p) p (N - n)/N,   ATI(p) = n + (1 - Pa(p)) (N - n).}
#' The **average outgoing quality limit** (AOQL) is the maximum of AOQ over `p`.
#'
#' @param plan A [sampling_plan()] with a lot size `N`.
#' @param p Fractions defective of the incoming lots.
#' @return `aoq()` and `ati()` return numeric vectors the same length as `p`. `aoql()` returns a
#'   named numeric vector with the incoming fraction defective `p` at which the AOQ is largest and
#'   the AOQL itself (`aoq`).
#' @examples
#' plan <- sampling_plan(50, 1, N = 5000)
#' aoq(plan, c(0.01, 0.02, 0.05))
#' ati(plan, c(0.01, 0.02, 0.05))
#' aoql(plan)
#' @name rectifying
NULL

#' @rdname rectifying
#' @export
aoq <- function(plan, p) {
  check_rectifying(plan)
  plan_pa(plan, p) * p * (plan$N - plan$n) / plan$N
}

#' @rdname rectifying
#' @export
ati <- function(plan, p) {
  check_rectifying(plan)
  plan$n + (1 - plan_pa(plan, p)) * (plan$N - plan$n)
}

#' @rdname rectifying
#' @export
aoql <- function(plan) {
  check_rectifying(plan)
  grid <- seq(0, 0.5, by = 0.0005)
  v <- aoq(plan, grid)
  c(p = grid[which.max(v)], aoq = max(v))
}

#' Find the smallest single sampling plan that meets two risks
#'
#' Searches for the single sampling plan `(n, c)` with the smallest sample size such that a lot
#' with fraction defective `aql` is accepted with probability at least `1 - alpha` and a lot with
#' fraction defective `ltpd` is accepted with probability at most `beta`.
#'
#' @param aql Acceptable quality level, as a fraction defective.
#' @param alpha Producer's risk: the largest acceptable probability of rejecting a lot at the
#'   AQL. Default 0.05.
#' @param ltpd Lot tolerance percent defective, as a fraction (greater than `aql`).
#' @param beta Consumer's risk: the largest acceptable probability of accepting a lot at the
#'   LTPD. Default 0.10.
#' @param N Lot size, or `NULL` (the default) for the binomial approximation.
#' @param n_max Largest sample size to consider. Default 2000.
#' @return A one-row tibble with the plan (`n`, `c`) and its `pa_aql`, `pa_ltpd`,
#'   `producers_risk`, and `consumers_risk`. An error if no plan with `n <= n_max` exists.
#' @examples
#' find_plan(aql = 0.01, ltpd = 0.05)
#' find_plan(aql = 0.01, ltpd = 0.05, N = 200)
#' @seealso [sampling_plan()], [plan_risks()]
#' @export
find_plan <- function(aql, alpha = 0.05, ltpd, beta = 0.10, N = NULL, n_max = 2000) {
  for (nm in c("aql", "alpha", "ltpd", "beta")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1 || is.na(v) || v <= 0 || v >= 1) {
      stop("`", nm, "` must be a single number between 0 and 1.", call. = FALSE)
    }
  }
  if (aql >= ltpd) stop("`aql` must be smaller than `ltpd`.", call. = FALSE)
  check_whole(n_max, "n_max", 2)
  if (!is.null(N)) n_max <- min(n_max, N)
  pa_at <- function(n, c, p) {
    if (is.null(N)) stats::pbinom(c, n, p) else {
      d <- round(N * p)
      stats::phyper(c, d, N - d, n)
    }
  }
  for (c in 0:(n_max - 1)) {
    n_vec <- (c + 1):n_max
    ok_ltpd <- which(pa_at(n_vec, c, ltpd) <= beta)
    if (length(ok_ltpd) == 0) next
    n <- n_vec[ok_ltpd[1]]
    if (pa_at(n, c, aql) >= 1 - alpha) {
      return(plan_risks(sampling_plan(n, c, N), aql, ltpd) |>
               (\(r) tibble::tibble(n = n, c = c, pa_aql = r$pa_aql, pa_ltpd = r$pa_ltpd,
                                    producers_risk = r$producers_risk, consumers_risk = r$consumers_risk))())
    }
  }
  stop("No plan with n <= ", n_max, " meets both risks.", call. = FALSE)
}

#' Simulate the acceptance probability of a sampling plan
#'
#' Simulates inspecting many lots with a given fraction defective and returns the proportion of
#' lots accepted, which estimates the OC curve. The sample is drawn from a lot of size `N`
#' (hypergeometric) or from an infinitely large lot (binomial).
#'
#' @param plan A [sampling_plan()].
#' @param p Fractions defective at which to simulate.
#' @param reps Number of simulated lots at each `p`. Default 10000.
#' @return A tibble with columns `p`, `pa_sim` (the proportion of lots accepted), and `se` (the
#'   standard error of that proportion).
#' @examples
#' set.seed(1)
#' sim_acceptance(sampling_plan(50, 1, N = 5000), p = c(0.01, 0.05), reps = 2000)
#' @seealso [oc_curve()]
#' @export
sim_acceptance <- function(plan, p, reps = 10000) {
  check_plan(plan)
  check_whole(reps, "reps", 1)
  if (!is.numeric(p) || anyNA(p) || any(p < 0) || any(p > 1)) stop("`p` must be fractions between 0 and 1.", call. = FALSE)
  pa_sim <- vapply(p, function(pi) {
    defects <- if (is.null(plan$N)) stats::rbinom(reps, plan$n, pi) else {
      d <- round(plan$N * pi)
      stats::rhyper(reps, d, plan$N - d, plan$n)
    }
    mean(defects <= plan$c)
  }, numeric(1))
  tibble::tibble(p = p, pa_sim = pa_sim, se = sqrt(pa_sim * (1 - pa_sim) / reps))
}
