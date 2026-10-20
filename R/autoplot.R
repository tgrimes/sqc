#' Draw a control chart with ggplot2
#'
#' Plots a `qc_chart` object: the plotted statistic in time order, the center line (solid),
#' the control limits (dashed), and signals in color (red for a point beyond the limits, and
#' orange for a one-sided run when `run_length` was set in [qc_chart()]).
#'
#' @param object A [qc_chart()] object.
#' @param title,subtitle Plot title and subtitle (character strings or plotmath expressions).
#'   The default title names the chart.
#' @param xlab,ylab Axis labels.
#' @param digits Number of decimals in the UCL, CL, and LCL labels.
#' @param ... Not used.
#' @return A ggplot object, so it can be modified with the usual `+` syntax.
#' @examples
#' q <- qc_chart(filler, type = "xbar", value = fill_ml, sample = cycle)
#' autoplot(q, subtitle = "Four-head filler, one bottle from each head")
#' @export
autoplot.qc_chart <- function(object, title = NULL, subtitle = NULL,
                              xlab = "Sample number (time order)", ylab = NULL, digits = 2, ...) {
  if (is.null(title)) {
    title <- if (object$type == "xbar") expression(bar(X) * " chart") else "R chart"
  }
  if (is.null(ylab)) ylab <- if (object$type == "xbar") "Sample mean" else "Sample range"
  has_runs <- is.finite(object$run_length)
  run_label <- paste0("Run of ", object$run_length, " on one side")
  d <- as.data.frame(object)
  d$status[d$status == "Run"] <- run_label
  lims <- unname(object$limits)
  m <- object$m
  pad <- 0.08 * diff(range(c(lims, d$statistic)))
  ylim <- c(min(lims, d$statistic) - pad, max(lims, d$statistic) + pad)
  edge <- data.frame(y = c(lims[2], object$center, lims[1]),
                     label = sprintf(paste0(c("UCL", "CL", "LCL"), " = %.", digits, "f"),
                                     c(lims[2], object$center, lims[1])))
  by <- max(1, 5 * ceiling(m / 50))
  breaks <- unique(c(1, if (m >= by) seq(by, m, by = by)))
  colors <- c("No signal" = "black", "Beyond limits" = "red")
  if (has_runs) colors[run_label] <- "orange"

  ggplot2::ggplot(d, ggplot2::aes(x = .data$sample, y = .data$statistic)) +
    ggplot2::geom_hline(yintercept = object$center) +
    ggplot2::geom_hline(yintercept = lims, linetype = "dashed") +
    ggplot2::geom_line(color = "grey60") +
    ggplot2::geom_point(ggplot2::aes(color = .data$status), size = 2.2) +
    ggplot2::geom_text(data = edge, ggplot2::aes(x = m + 0.7, y = .data$y, label = .data$label),
                       hjust = 0, size = 3.6, inherit.aes = FALSE) +
    ggplot2::scale_color_manual(values = colors, breaks = c("Beyond limits", if (has_runs) run_label), name = NULL) +
    ggplot2::scale_x_continuous(breaks = breaks) +
    ggplot2::coord_cartesian(xlim = c(0.5, m + 0.5), ylim = ylim, expand = FALSE, clip = "off") +
    ggplot2::labs(x = xlab, y = ylab, title = title, subtitle = subtitle) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.margin = ggplot2::margin(6, 70, 6, 6), legend.position = "bottom")
}

#' @export
plot.qc_chart <- function(x, ...) {
  print(autoplot.qc_chart(x, ...))
  invisible(x)
}

# Re-export the generic so that `autoplot(q)` works after `library(sqc)` alone.
#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot
