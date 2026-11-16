#' Gage R&R: a crossed two-factor random effects study (ANOVA estimators)
#'
#' Analyzes a **crossed, balanced** measurement system study in which every part (sample, bag, patient)
#' is measured by every operator (technician, analyst), with `r` replicate measurements for each
#' part and operator. The model is the two-factor random effects model with interaction,
#' \deqn{y_{ijk} = \mu + P_i + O_j + (PO)_{ij} + e_{ijk},}{y_ijk = mu + P_i + O_j + (PO)_ij + e_ijk,}
#' with independent normal effects with variances \eqn{\sigma_P^2}{sigma_P^2} (part-to-part),
#' \eqn{\sigma_O^2}{sigma_O^2} (operator), \eqn{\sigma_{PO}^2}{sigma_PO^2} (part by operator
#' interaction), and \eqn{\sigma_e^2}{sigma_e^2} (repeatability).
#'
#' The variance components are estimated by the **ANOVA (method of moments) estimators**
#' \deqn{\hat\sigma_e^2 = MS_E, \quad \hat\sigma_{PO}^2 = \frac{MS_{PO} - MS_E}{r}, \quad
#' \hat\sigma_O^2 = \frac{MS_O - MS_{PO}}{pr}, \quad \hat\sigma_P^2 = \frac{MS_P - MS_{PO}}{or},}{sigma_e^2 = MSE, sigma_PO^2 = (MSPO - MSE)/r, sigma_O^2 = (MSO - MSPO)/(p r), sigma_P^2 = (MSP - MSPO)/(o r),}
#' where \eqn{p}, \eqn{o}, \eqn{r} are the numbers of parts, operators, and replicates. Parts and operators are
#' tested against the interaction mean square, and the interaction against the error mean square.
#'
#' With `interaction = "auto"` (the AIAG convention) the interaction is dropped when its p-value is
#' above `alpha` (0.25): the interaction and error sums of squares are pooled, the error mean square
#' becomes \eqn{(SS_{PO} + SS_E)/(df_{PO} + df_E)}{(SSPO + SSE)/(dfPO + dfE)}, and
#' \eqn{\hat\sigma_{PO}^2}{sigma_PO^2} is 0. Use `"yes"` or `"no"` to force the choice.
#'
#' A negative estimate of a component is shown as it is in the `variance` column and treated as 0 in the
#' standard deviations, the percentages, and the totals.
#'
#' The Gage R&R quantities are
#' * **repeatability** \eqn{\sigma_e}{sigma_e};
#' * **reproducibility** \eqn{\sqrt{\sigma_O^2 + \sigma_{PO}^2}}{sqrt(sigma_O^2 + sigma_PO^2)};
#' * **gage R&R** \eqn{\sigma_{GRR}^2 = \sigma_e^2 + \sigma_O^2 + \sigma_{PO}^2}{sigma_GRR^2 = sigma_e^2 + sigma_O^2 + sigma_PO^2};
#' * **total variation** \eqn{\sigma_T^2 = \sigma_{GRR}^2 + \sigma_P^2}{sigma_T^2 = sigma_GRR^2 + sigma_P^2};
#' * the **percent contribution** \eqn{100\,\sigma_{GRR}^2/\sigma_T^2}{100 sigma_GRR^2/sigma_T^2};
#' * the **percent study variation** (the "%GRR") \eqn{100\,\sigma_{GRR}/\sigma_T}{100 sigma_GRR/sigma_T};
#' * the **percent of tolerance** \eqn{100 \cdot 6\sigma_{GRR}/\text{tolerance}}{100 * 6 sigma_GRR/tolerance} when `tolerance` is given;
#' * the **number of distinct categories** \eqn{\lfloor\sqrt{2}\,\sigma_P/\sigma_{GRR}\rfloor}{floor(sqrt(2) sigma_P/sigma_GRR)}.
#'
#' The usual (AIAG) guide for the %GRR is below 10 percent acceptable, 10 to 30 percent marginal, and
#' above 30 percent unacceptable, with at least 5 distinct categories.
#'
#' @param data A data frame with one row per measurement.
#' @param value,part,operator The columns holding the measurement, the part, and the operator
#'   (bare names, as in dplyr). Every part-by-operator cell must have the same number \eqn{r \ge 2}{r >= 2}
#'   of replicates, with at least 2 parts and 2 operators.
#' @param interaction `"auto"` (default), `"yes"` (always keep the interaction), or `"no"` (always pool it).
#' @param alpha The significance level for dropping the interaction when `interaction = "auto"`. Default 0.25.
#' @param tolerance Optional specification width (upper minus lower limit) for the percent of tolerance.
#' @return An object of class `gage_rr`: a list with `anova` (a tibble with the source, degrees of
#'   freedom, sums of squares, mean squares, F and p-value, using the model that was kept), `components` (a tibble
#'   with the variance, standard deviation, percent contribution, and percent study variation of each
#'   source), `summary` (a named list: `pct_grr`, `pct_contribution`, `ndc`, `sd_repeatability`,
#'   `sd_reproducibility`, `sd_grr`, `sd_part`, `sd_total`, `pct_tolerance`), `p_interaction` (the p-value of the
#'   interaction test), `interaction_used`, `negative` (names of the components with a negative estimate),
#'   `design` (`parts`, `operators`, `replicates`), and `alpha`.
#' @seealso [sim_gage()] to simulate the sampling distribution of the estimates, and [vc_oneway()] for
#'   the one-way model.
#' @examples
#' gage_rr(assay_study, value = hba1c, part = sample, operator = technician)
#'
#' # Force the interaction out of the model, and add a specification width
#' gage_rr(assay_study, value = hba1c, part = sample, operator = technician,
#'         interaction = "no", tolerance = 4)
#' @export
gage_rr <- function(data, value, part, operator, interaction = c("auto", "yes", "no"), alpha = 0.25,
                    tolerance = NULL) {
  interaction <- match.arg(interaction)
  check_probability(alpha, "alpha")
  if (!is.null(tolerance)) check_number(tolerance, "tolerance", positive = TRUE)
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  if (missing(value) || missing(part) || missing(operator)) {
    stop("Supply `value`, `part`, and `operator` (the measurement, part, and operator columns).", call. = FALSE)
  }
  y <- rlang::eval_tidy(rlang::enquo(value), data)
  pa <- rlang::eval_tidy(rlang::enquo(part), data)
  op <- rlang::eval_tidy(rlang::enquo(operator), data)
  if (!is.numeric(y) || anyNA(y) || anyNA(pa) || anyNA(op)) {
    stop("Measurements must be numeric, and measurements, parts, and operators must have no missing values.", call. = FALSE)
  }
  pa <- factor(pa, levels = unique(pa))
  op <- factor(op, levels = unique(op))
  p <- nlevels(pa); o <- nlevels(op)
  if (p < 2 || o < 2) stop("At least 2 parts and 2 operators are needed.", call. = FALSE)
  cell <- table(pa, op)
  if (length(unique(as.vector(cell))) != 1 || cell[1, 1] < 2) {
    stop("The design must be balanced with the same number of replicates (at least 2) for every part and operator.", call. = FALSE)
  }
  r <- cell[[1, 1]]

  cell_mean <- tapply(y, list(pa, op), mean)
  part_mean <- rowMeans(cell_mean)
  op_mean <- colMeans(cell_mean)
  grand <- mean(y)
  ss <- c(part = o * r * sum((part_mean - grand)^2),
          operator = p * r * sum((op_mean - grand)^2),
          interaction = r * sum((cell_mean - outer(part_mean, rep(1, o)) - outer(rep(1, p), op_mean) + grand)^2),
          error = sum((y - cell_mean[cbind(as.integer(pa), as.integer(op))])^2))
  gage_from_ss(ss, p, o, r, interaction, alpha, tolerance)
}

# Everything after the sums of squares: the ANOVA table, the components, and the Gage R&R summary.
gage_from_ss <- function(ss, p, o, r, interaction, alpha, tolerance) {
  df <- c(part = p - 1, operator = o - 1, interaction = (p - 1) * (o - 1), error = p * o * (r - 1))
  ms <- ss / df
  p_int <- stats::pf(ms[["interaction"]] / ms[["error"]], df[["interaction"]], df[["error"]], lower.tail = FALSE)
  keep <- switch(interaction, yes = TRUE, no = FALSE, auto = p_int <= alpha)
  est <- gage_estimates(ms, ss, df, p, o, r, keep)

  denom <- if (keep) ms[["interaction"]] else (ss[["interaction"]] + ss[["error"]]) / (df[["interaction"]] + df[["error"]])
  denom_df <- if (keep) df[["interaction"]] else df[["interaction"]] + df[["error"]]
  f_part <- ms[["part"]] / denom; f_op <- ms[["operator"]] / denom
  f_int <- ms[["interaction"]] / ms[["error"]]
  p_part <- stats::pf(f_part, df[["part"]], denom_df, lower.tail = FALSE)
  p_op <- stats::pf(f_op, df[["operator"]], denom_df, lower.tail = FALSE)
  anova <- tibble::tibble(
    source = c("Part", "Operator", "Part x operator", "Repeatability (error)"),
    df = unname(df), ss = unname(ss), ms = unname(ms),
    f = c(f_part, f_op, f_int, NA_real_),
    p = c(p_part, p_op, p_int, NA_real_)
  )

  v <- est
  v[v < 0] <- 0                                  # negative estimates are treated as 0
  grr <- v[["repeatability"]] + v[["operator"]] + v[["interaction"]]
  total <- grr + v[["part"]]
  repro <- v[["operator"]] + v[["interaction"]]
  var_tab <- c(est[["repeatability"]], repro, est[["operator"]], est[["interaction"]], grr, est[["part"]], total)
  var_used <- c(v[["repeatability"]], repro, v[["operator"]], v[["interaction"]], grr, v[["part"]], total)
  names_tab <- c("Repeatability", "Reproducibility", "  Operator", "  Part x operator", "Gage R&R", "Part-to-part", "Total variation")
  components <- tibble::tibble(
    component = names_tab, variance = var_tab, sd = sqrt(var_used),
    percent_contribution = 100 * var_used / total, percent_study = 100 * sqrt(var_used) / sqrt(total)
  )
  summary <- list(
    pct_grr = 100 * sqrt(grr / total), pct_contribution = 100 * grr / total,
    ndc = floor(sqrt(2) * sqrt(v[["part"]] / grr)),
    sd_repeatability = sqrt(v[["repeatability"]]), sd_reproducibility = sqrt(repro), sd_grr = sqrt(grr),
    sd_part = sqrt(v[["part"]]), sd_total = sqrt(total),
    pct_tolerance = if (is.null(tolerance)) NA_real_ else 100 * 6 * sqrt(grr) / tolerance
  )
  structure(
    list(anova = anova, components = components, summary = summary, p_interaction = p_int, interaction_used = keep,
         negative = names(est)[est < 0], design = list(parts = p, operators = o, replicates = r), alpha = alpha),
    class = "gage_rr"
  )
}

# ANOVA estimates of the four variance components, with the interaction kept or pooled with the error.
gage_estimates <- function(ms, ss, df, p, o, r, keep) {
  if (keep) {
    c(repeatability = ms[["error"]],
      operator = (ms[["operator"]] - ms[["interaction"]]) / (p * r),
      interaction = (ms[["interaction"]] - ms[["error"]]) / r,
      part = (ms[["part"]] - ms[["interaction"]]) / (o * r))
  } else {
    pooled <- (ss[["interaction"]] + ss[["error"]]) / (df[["interaction"]] + df[["error"]])
    c(repeatability = pooled,
      operator = (ms[["operator"]] - pooled) / (p * r),
      interaction = 0,
      part = (ms[["part"]] - pooled) / (o * r))
  }
}

#' @export
print.gage_rr <- function(x, digits = 4, ...) {
  d <- x$design
  cat("Gage R&R (crossed, ANOVA estimators): ", d$parts, " parts, ", d$operators, " operators, ", d$replicates,
      " replicates\n\n", sep = "")
  an <- as.data.frame(x$anova)
  print(format_table(an, digits), row.names = FALSE)
  cat("\nPart x operator interaction: p = ", format(signif(x$p_interaction, 3)),
      if (x$interaction_used) " (kept in the model)\n\n" else " (pooled with the error)\n\n", sep = "")
  comp <- as.data.frame(x$components)
  comp$percent_contribution <- round(comp$percent_contribution, 1)
  comp$percent_study <- round(comp$percent_study, 1)
  print(format_table(comp, digits), row.names = FALSE)
  if (length(x$negative) > 0) {
    cat("\nNegative ANOVA estimate for: ", paste(x$negative, collapse = ", "),
        "; the sd, percentages, and totals use 0.\n", sep = "")
  }
  s <- x$summary
  rating <- if (s$pct_grr < 10) "acceptable" else if (s$pct_grr <= 30) "marginal" else "unacceptable"
  cat("\n%GRR (study variation) = ", format(round(s$pct_grr, 1)), "  (", rating, " by the 10/30 guide)",
      "\nNumber of distinct categories = ", s$ndc, "\n", sep = "")
  if (!is.na(s$pct_tolerance)) cat("% of tolerance = ", format(round(s$pct_tolerance, 1)), "\n", sep = "")
  invisible(x)
}

#' Simulate the Gage R&R estimates of a crossed two-factor study
#'
#' Generates `nsim` balanced crossed studies from the two-factor random effects model of [gage_rr()]
#' with `p` parts, `o` operators, and `r` replicates, applies the ANOVA estimators to each, and returns the
#' estimated components, the %GRR, and the number of distinct categories of each study. The spread
#' of the %GRR over the simulated studies shows how precisely a study of this size measures the
#' gage, and how often it would put the gage in the wrong category; the share of negative estimates of the
#' operator component shows what happens with few operators.
#'
#' Use [set.seed()] before the call for a reproducible result.
#'
#' @param p,o,r Number of parts, operators, and replicates (each at least 2).
#' @param sigma_p,sigma_o,sigma_po,sigma_e Standard deviations of the part effects, operator effects,
#'   interaction effects (all non-negative), and repeatability error (positive).
#' @param nsim Number of simulated studies. Default 10000.
#' @param interaction,alpha As in [gage_rr()]. The default `"auto"` applies the 0.25 rule to each study.
#' @return A tibble with one row per study and columns `sigma2_p`, `sigma2_o`, `sigma2_po`, `sigma2_e` (the ANOVA
#'   estimates, which can be negative), `p_interaction`, `sd_grr`, `pct_grr`, and `ndc`.
#' @seealso [gage_rr()]
#' @examples
#' set.seed(1)
#' s <- sim_gage(p = 10, o = 3, r = 2, sigma_p = 1.4, sigma_o = 0.2, sigma_po = 0.08,
#'               sigma_e = 0.15, nsim = 2000)
#' quantile(s$pct_grr, c(0.05, 0.5, 0.95))   # true %GRR is about 18
#' mean(s$sigma2_o < 0)                      # share of negative operator estimates
#' @export
sim_gage <- function(p, o, r, sigma_p, sigma_o, sigma_po, sigma_e, nsim = 10000,
                     interaction = c("auto", "yes", "no"), alpha = 0.25) {
  interaction <- match.arg(interaction)
  check_count(p, "p", minimum = 2); check_count(o, "o", minimum = 2); check_count(r, "r", minimum = 2)
  check_count(nsim, "nsim", minimum = 1)
  check_probability(alpha, "alpha")
  for (nm in c("sigma_p", "sigma_o", "sigma_po")) {
    v <- get(nm)
    if (!(is.numeric(v) && length(v) == 1 && is.finite(v) && v >= 0)) {
      stop("`", nm, "` must be a single non-negative finite number.", call. = FALSE)
    }
  }
  check_number(sigma_e, "sigma_e", positive = TRUE)

  P <- matrix(stats::rnorm(p * nsim), p, nsim) * sigma_p
  O <- matrix(stats::rnorm(o * nsim), o, nsim) * sigma_o
  PO <- array(stats::rnorm(p * o * nsim), c(p, o, nsim)) * sigma_po
  E <- array(stats::rnorm(r * p * o * nsim), c(r, p, o, nsim)) * sigma_e
  P3 <- aperm(array(P, c(p, nsim, o)), c(1, 3, 2))
  O3 <- aperm(array(O, c(o, nsim, p)), c(3, 1, 2))
  Y <- E + rep(as.vector(P3 + O3 + PO), each = r)                          # r x p x o x nsim
  cell <- colMeans(Y)                                                       # p x o x nsim
  ss_e <- colSums(matrix(colSums((Y - rep(as.vector(cell), each = r))^2), p * o, nsim))
  pm <- colMeans(aperm(cell, c(2, 1, 3)))                                   # part means, p x nsim
  om <- colMeans(cell)                                                       # operator means, o x nsim
  gm <- colMeans(pm)
  ss_p <- o * r * colSums((pm - rep(gm, each = p))^2)
  ss_o <- p * r * colSums((om - rep(gm, each = o))^2)
  dev <- cell - aperm(array(pm, c(p, nsim, o)), c(1, 3, 2)) - aperm(array(om, c(o, nsim, p)), c(3, 1, 2)) +
    rep(gm, each = p * o)
  ss_po <- r * colSums(matrix(dev^2, p * o, nsim))

  df <- c(part = p - 1, operator = o - 1, interaction = (p - 1) * (o - 1), error = p * o * (r - 1))
  ms_p <- ss_p / df[["part"]]; ms_o <- ss_o / df[["operator"]]
  ms_po <- ss_po / df[["interaction"]]; ms_e <- ss_e / df[["error"]]
  p_int <- stats::pf(ms_po / ms_e, df[["interaction"]], df[["error"]], lower.tail = FALSE)
  keep <- switch(interaction, yes = rep(TRUE, nsim), no = rep(FALSE, nsim), auto = p_int <= alpha)
  pooled <- (ss_po + ss_e) / (df[["interaction"]] + df[["error"]])
  s2_e <- ifelse(keep, ms_e, pooled)
  s2_po <- ifelse(keep, (ms_po - ms_e) / r, 0)
  den <- ifelse(keep, ms_po, pooled)
  s2_o <- (ms_o - den) / (p * r)
  s2_p <- (ms_p - den) / (o * r)
  grr <- pmax(0, s2_e) + pmax(0, s2_o) + pmax(0, s2_po)
  total <- grr + pmax(0, s2_p)
  tibble::tibble(sigma2_p = s2_p, sigma2_o = s2_o, sigma2_po = s2_po, sigma2_e = s2_e, p_interaction = p_int,
                 sd_grr = sqrt(grr), pct_grr = 100 * sqrt(grr / total),
                 ndc = floor(sqrt(2) * sqrt(pmax(0, s2_p) / grr)))
}

check_probability <- function(x, name) {
  if (!(is.numeric(x) && length(x) == 1 && is.finite(x) && x > 0 && x < 1)) {
    stop("`", name, "` must be a single number between 0 and 1.", call. = FALSE)
  }
  invisible(x)
}
