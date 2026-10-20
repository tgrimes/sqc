#' sqc: companion tools for a Statistical Quality Control course
#'
#' The package provides three things:
#' * control chart constants: [d2()], [d3()], [c4()], and the table [qc_constants()];
#' * Shewhart charts: [qc_chart()] computes the center line, limits, and signals, and
#'   [autoplot.qc_chart()] draws the chart with ggplot2;
#' * example datasets: [filler] and [filler_rotating].
#'
#' @keywords internal
#' @importFrom rlang .data
"_PACKAGE"
