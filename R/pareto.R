#' Pareto chart
#'
#' Draws a **Pareto chart**: the categories sorted by count (bars) with the cumulative percentage (line), to show which few
#' categories account for most of the problems.
#'
#' @param counts A named numeric vector of counts, or a data frame with the columns `category` and `count`.
#' @param top Keep only the `top` largest categories and combine the rest in one bar "Other" (placed last). Default `NULL` keeps all.
#' @param title,subtitle Plot title and subtitle.
#' @param xlab,ylab Axis labels.
#' @return A ggplot object.
#' @examples
#' pareto_chart(c(CLABSI = 41, CAUTI = 78, `SSI colon` = 23, `SSI hysterectomy` = 9, MRSA = 12,
#'                `C. difficile` = 64))
#' @export
pareto_chart <- function(counts, top = NULL, title = "Pareto chart", subtitle = NULL, xlab = NULL, ylab = "Count") {
  if (is.data.frame(counts)) {
    if (!all(c("category", "count") %in% names(counts))) stop("A data frame needs the columns `category` and `count`.", call. = FALSE)
    counts <- stats::setNames(counts$count, counts$category)
  }
  if (!is.numeric(counts) || is.null(names(counts)) || anyNA(counts) || any(counts < 0) || length(counts) < 2) {
    stop("`counts` must be a named numeric vector of at least 2 non-negative counts.", call. = FALSE)
  }
  counts <- sort(counts, decreasing = TRUE)
  if (!is.null(top)) {
    check_count(top, "top", minimum = 1)
    if (top < length(counts)) counts <- c(counts[seq_len(top)], Other = sum(counts[-seq_len(top)]))
  }
  total <- sum(counts)
  d <- tibble::tibble(category = factor(names(counts), levels = names(counts)), count = as.numeric(counts),
                      cum = cumsum(as.numeric(counts)) / total)
  ggplot2::ggplot(d, ggplot2::aes(x = .data$category)) +
    ggplot2::geom_col(ggplot2::aes(y = .data$count), fill = "#4c78a8", width = 0.75) +
    ggplot2::geom_line(ggplot2::aes(y = .data$cum * total, group = 1), color = "#e45756", linewidth = 0.9) +
    ggplot2::geom_point(ggplot2::aes(y = .data$cum * total), color = "#e45756", size = 2) +
    ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(~ . / total, name = "Cumulative percentage", labels = function(v) paste0(round(100 * v), "%"))) +
    ggplot2::labs(x = xlab, y = ylab, title = title, subtitle = subtitle) +
    ggplot2::theme_bw() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}
