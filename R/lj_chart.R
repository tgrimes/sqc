#' Levey-Jennings chart for control results
#'
#' Computes a Levey-Jennings chart: the results of one or more control materials (levels), run by
#' run, compared with the established mean and standard deviation of each material, and checked
#' with the Westgard rules. A Levey-Jennings chart is an individuals chart with known parameters
#' and lines at the mean and at 1, 2, and 3 SD, the standard tool of clinical laboratory quality
#' control. Draw it with [autoplot.lj_chart()].
#'
#' Unlike [qc_chart()], which estimates the center line and limits from subgroup data, the mean
#' and SD here are given: they come from earlier data on the control material (for example, 20 or
#' more runs collected when the material was first put into use).
#'
#' @param x A numeric matrix or data frame of control results with one row per run (in time order)
#'   and one column per control level, or a numeric vector for a single level.
#' @param mean,sd The established mean and standard deviation of each control level: a single
#'   number (used for every level) or one number per level.
#' @param rules The rules used to reject a run: any of `"1-2s"`, `"1-3s"`, `"2-2s"`, `"R-4s"`,
#'   `"4-1s"`, `"10x"`, or `"multirule"` (the default). See [apply_rules()].
#' @param level_labels Optional labels for the levels (default `"Level 1"`, `"Level 2"`, ...).
#' @return An object of class `lj_chart`, a list with the control results (`values`), their
#'   z-scores (`z`), the established `mean` and `sd`, the `rules`, the rule `violations` (see
#'   [apply_rules()]), the indices of the `rejected` runs, and `level_labels`.
#' @examples
#' set.seed(1)
#' x <- cbind(rnorm(30, 100, 3), rnorm(30, 250, 6))
#' x[19:30, ] <- x[19:30, ] + rep(c(4.5, 9), each = 12)   # a 1.5 SD shift from run 19
#' q <- lj_chart(x, mean = c(100, 250), sd = c(3, 6))
#' q
#' @seealso [apply_rules()], [sim_run_length()]
#' @export
lj_chart <- function(x, mean, sd, rules = "multirule", level_labels = NULL) {
  x <- as.matrix(x)
  if (!is.numeric(x) || anyNA(x)) stop("`x` must be numeric with no missing values.", call. = FALSE)
  n_levels <- ncol(x)
  for (nm in c("mean", "sd")) {
    v <- get(nm)
    if (!is.numeric(v) || anyNA(v) || !(length(v) %in% c(1, n_levels))) {
      stop("`", nm, "` must be a number, or one number per control level.", call. = FALSE)
    }
  }
  if (any(sd <= 0)) stop("`sd` must be positive.", call. = FALSE)
  mean <- rep_len(mean, n_levels)
  sd <- rep_len(sd, n_levels)
  if (is.null(level_labels)) level_labels <- paste("Level", seq_len(n_levels))
  if (length(level_labels) != n_levels) stop("`level_labels` must have one label per level.", call. = FALSE)
  rules <- expand_rules(rules)
  z <- sweep(sweep(x, 2, mean), 2, sd, "/")
  violations <- apply_rules(z, rules)
  structure(list(values = x, z = z, mean = mean, sd = sd, rules = rules, violations = violations,
                 rejected = sort(unique(violations$run)), level_labels = level_labels,
                 m = nrow(x), n_levels = n_levels),
            class = "lj_chart")
}

#' @export
print.lj_chart <- function(x, ...) {
  cat("<lj_chart> ", x$m, " runs, ", x$n_levels, " control level(s)\n", sep = "")
  cat("  rules:         ", paste(x$rules, collapse = ", "), "\n", sep = "")
  cat("  runs rejected: ", length(x$rejected), if (length(x$rejected) > 0)
        paste0(" (first at run ", min(x$rejected), ")") else "", "\n", sep = "")
  if (nrow(x$violations) > 0) {
    counts <- table(factor(x$violations$rule, levels = x$rules))
    counts <- counts[counts > 0]
    cat("  violations:    ", paste0(names(counts), ": ", counts, collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}

#' @export
as.data.frame.lj_chart <- function(x, ...) {
  data.frame(run = rep(seq_len(x$m), times = x$n_levels),
             level = rep(seq_len(x$n_levels), each = x$m),
             value = as.vector(x$values), z = as.vector(x$z),
             rejected = rep(seq_len(x$m) %in% x$rejected, times = x$n_levels))
}

#' Draw a Levey-Jennings chart with ggplot2
#'
#' One panel per control level. Each panel shows the results in run order with lines at the
#' mean (solid) and at plus and minus 1 SD (dotted), 2 SD (dashed amber), and 3 SD (dashed red).
#' A rejected run is shaded in all panels and labeled with the rule(s) that fired.
#'
#' @param object An [lj_chart()] object.
#' @param title,subtitle Plot title and subtitle. The default title is "Levey-Jennings chart" and
#'   the default subtitle lists the rules.
#' @param xlab,ylab Axis labels.
#' @param show_rules If `TRUE` (default), label each rejected run with the rules that fired.
#' @param ... Not used.
#' @return A ggplot object.
#' @examples
#' set.seed(1)
#' x <- cbind(rnorm(30, 100, 3), rnorm(30, 250, 6))
#' autoplot(lj_chart(x, mean = c(100, 250), sd = c(3, 6)))
#' @export
autoplot.lj_chart <- function(object, title = "Levey-Jennings chart", subtitle = NULL,
                              xlab = "Run", ylab = "Control result", show_rules = TRUE, ...) {
  lv <- factor(object$level_labels, levels = object$level_labels)
  m <- object$m
  pts <- data.frame(run = rep(seq_len(m), times = object$n_levels),
                    level = rep(lv, each = m), value = as.vector(object$values))
  ks <- -3:3
  kinds <- c("Mean", "\u00b11 SD", "\u00b12 SD", "\u00b13 SD")
  lines <- do.call(rbind, lapply(seq_len(object$n_levels), function(l) {
    data.frame(level = lv[l], k = ks, y = object$mean[l] + ks * object$sd[l],
               kind = factor(kinds[abs(ks) + 1], levels = kinds),
               label = ifelse(ks == 0, "Mean", sprintf("%+d SD", ks)))
  }))
  range_pts <- do.call(rbind, lapply(seq_len(object$n_levels), function(l) {
    data.frame(level = lv[l], run = 1, value = object$mean[l] + c(-3.5, 3.5) * object$sd[l])
  }))
  if (is.null(subtitle)) subtitle <- paste("Rules:", paste(object$rules, collapse = ", "))
  by <- max(1, 5 * ceiling(m / 50))
  breaks <- unique(c(1, if (m >= by) seq(by, m, by = by)))

  p <- ggplot2::ggplot(pts, ggplot2::aes(x = .data$run, y = .data$value)) +
    ggplot2::geom_blank(data = range_pts) +
    ggplot2::geom_hline(data = lines, ggplot2::aes(yintercept = .data$y, color = .data$kind, linetype = .data$kind)) +
    ggplot2::scale_color_manual(values = c(Mean = "black", "\u00b11 SD" = "grey60",
                                           "\u00b12 SD" = "#d97706", "\u00b13 SD" = "red3"), name = NULL) +
    ggplot2::scale_linetype_manual(values = c(Mean = "solid", "\u00b11 SD" = "dotted",
                                              "\u00b12 SD" = "dashed", "\u00b13 SD" = "dashed"), name = NULL)
  if (length(object$rejected) > 0) {
    bands <- do.call(rbind, lapply(seq_len(object$n_levels), function(l) {
      data.frame(level = lv[l], xmin = object$rejected - 0.5, xmax = object$rejected + 0.5)
    }))
    p <- p + ggplot2::geom_rect(data = bands, ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = -Inf, ymax = Inf),
                                fill = "red", alpha = 0.12, inherit.aes = FALSE)
    if (show_rules) {
      v <- object$violations
      txt <- vapply(object$rejected, function(r) paste(unique(v$rule[v$run == r]), collapse = ", "), character(1))
      labs_df <- data.frame(level = lv[1], run = object$rejected, label = txt)
      p <- p + ggplot2::geom_text(data = labs_df, ggplot2::aes(x = .data$run, y = Inf, label = .data$label),
                                  angle = 90, hjust = 1.05, vjust = -0.4, size = 2.8, color = "red3", inherit.aes = FALSE)
    }
  }
  p +
    ggplot2::geom_line(color = "grey55") +
    ggplot2::geom_point(size = 1.9) +
    ggplot2::geom_text(data = lines, ggplot2::aes(x = m + 1.2, y = .data$y, label = .data$label),
                       hjust = 0, size = 2.8, inherit.aes = FALSE) +
    ggplot2::facet_wrap(~ level, ncol = 1, scales = "free_y") +
    ggplot2::scale_x_continuous(breaks = breaks, expand = c(0, 0)) +
    ggplot2::coord_cartesian(xlim = c(0.5, m + 0.5), clip = "off") +
    ggplot2::labs(x = xlab, y = ylab, title = title, subtitle = subtitle) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.margin = ggplot2::margin(6, 55, 6, 6), legend.position = "bottom")
}
