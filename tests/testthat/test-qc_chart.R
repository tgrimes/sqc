x <- matrix(c(10, 12, 11, 13,
              11, 11, 12, 10,
              9, 10, 12, 11,
              12, 13, 12, 14), ncol = 4, byrow = TRUE)
rbar <- mean(c(3, 2, 3, 2))   # subgroup ranges
xbarbar <- mean(rowMeans(x))

test_that("xbar chart matches the hand calculation", {
  q <- qc_chart(x, type = "xbar")
  sigma_hat <- rbar / d2(4)
  expect_equal(q$center, xbarbar)
  expect_equal(q$sigma, sigma_hat)
  expect_equal(unname(q$limits), xbarbar + c(-1, 1) * 3 * sigma_hat / sqrt(4))
  expect_equal(unname(q$limits), xbarbar + c(-1, 1) * qc_constants(4)$A2 * rbar)
  expect_equal(q$statistics, rowMeans(x))
  expect_equal(c(q$m, q$n), c(4, 4))
})

test_that("R chart matches the hand calculation", {
  q <- qc_chart(x, type = "r")
  cons <- qc_constants(4)
  expect_equal(q$center, rbar)
  expect_equal(q$sigma, rbar / d2(4))
  expect_equal(unname(q$limits), c(cons$D3 * rbar, cons$D4 * rbar))
})

test_that("known center and sigma are used as given", {
  q <- qc_chart(x, type = "xbar", center = 11, sigma = 1)
  expect_equal(q$center, 11)
  expect_equal(unname(q$limits), 11 + c(-1, 1) * 1.5)
  r <- qc_chart(x, type = "r", sigma = 1)
  expect_equal(r$center, d2(4))
  expect_equal(r$limits[["lcl"]], max(0, d2(4) - 3 * d3(4)))
  r2 <- qc_chart(x, type = "r", center = 2)
  expect_equal(r2$sigma, 2 / d2(4))
})

test_that("a data frame gives the same chart as the matrix", {
  long <- data.frame(s = rep(1:4, each = 4), v = as.vector(t(x)))
  q_df <- qc_chart(long, type = "xbar", value = v, sample = s)
  q_mat <- qc_chart(x, type = "xbar")
  expect_equal(q_df$statistics, q_mat$statistics)
  expect_equal(q_df$limits, q_mat$limits)
  expect_equal(qc_chart(filler, "xbar", fill_ml, cycle)$m, 40)
})

test_that("signals: points beyond the limits and one-sided runs", {
  # six subgroups near the center, then one far above it
  near <- matrix(c(10, 10.5, 9.5, 10), nrow = 6, ncol = 4, byrow = TRUE)
  far <- c(30, 31, 30, 31)
  q <- qc_chart(rbind(near, far), type = "xbar", center = 10, sigma = 0.5)
  expect_equal(q$signals$beyond, 7)

  # nine subgroup means just above the center line form a run (but stay inside the limits);
  # the signal fires at the 7th point
  above <- matrix(10.3, nrow = 9, ncol = 4)
  q2 <- qc_chart(above, type = "xbar", center = 10, sigma = 1, run_length = 7)
  expect_equal(q2$signals$run, 7:9)   # flagged from the 7th point of the run
  expect_length(q2$signals$beyond, 0)

  # a run shorter than run_length is not flagged
  q3 <- qc_chart(above[1:6, ], type = "xbar", center = 10, sigma = 1, run_length = 7)
  expect_length(q3$signals$run, 0)
})

test_that("bad input gives clear errors", {
  expect_error(qc_chart(1:10), "matrix")
  expect_error(qc_chart(filler, "xbar"), "value.*sample")
  expect_error(qc_chart(x, "xbar", value = v), "only used with a data frame")
  expect_error(qc_chart(matrix(1:4, 4, 1)), "at least 2")
  ragged <- data.frame(s = c(1, 1, 1, 2, 2), v = 1:5)
  expect_error(qc_chart(ragged, "xbar", v, s), "same size")
  expect_error(qc_chart(x, sigma = -1), "positive")
})

test_that("print and autoplot work", {
  q <- qc_chart(x, type = "xbar")
  expect_output(print(q), "X-bar chart")
  p <- ggplot2::autoplot(q, title = "Test chart", subtitle = "sub")
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  expect_equal(p$labels$title, "Test chart")
  r <- ggplot2::autoplot(qc_chart(x, type = "r"))
  expect_equal(r$labels$title, "R chart")
})

test_that("the four-head filler shows the expected widening", {
  mixed <- qc_chart(subset(filler, cycle <= 20), "xbar", fill_ml, cycle)
  rational <- qc_chart(filler_rotating, "xbar", fill_ml, (cycle - 1) %/% 5 + 1)
  expect_gt(mixed$sigma, 1.5)       # head differences inflate sigma-hat
  expect_lt(rational$sigma, 1.3)    # within-head value is 1
  expect_gt(length(rational$signals$beyond), 0)
  expect_length(mixed$signals$beyond, 0)
})

test_that("subgroup ranges match a row-by-row calculation", {
  set.seed(5)
  m <- matrix(round(rnorm(60, sd = 3)), 12, 5)          # rounding creates ties; includes negatives
  expected <- apply(m, 1, function(r) max(r) - min(r))
  q <- qc_chart(m, type = "R")
  expect_equal(q$statistics, expected)
  expect_equal(q$center, mean(expected))
  expect_equal(qc_chart(m[, 1:2], type = "R")$statistics, abs(m[, 1] - m[, 2]))
})

test_that("invalid k is an error", {
  for (bad in list(-1, 0, NA_real_, Inf, c(2, 3), "3", NULL)) {
    expect_error(qc_chart(x, "xbar", k = bad), "`k` must be a single positive finite number")
  }
  expect_no_error(qc_chart(x, "xbar", k = 2.5))
  # limits move with k
  expect_equal(qc_chart(x, "xbar", k = 2)$limits[["ucl"]] - qc_chart(x, "xbar", k = 2)$center,
               2 * qc_chart(x, "xbar")$sigma / sqrt(4))
})

test_that("invalid run_length is an error, and Inf turns the run signal off", {
  for (bad in list(0, 1, 2.5, NA_real_, -Inf, "a", c(7, 8), NULL)) {
    expect_error(qc_chart(x, "xbar", run_length = bad), "`run_length` must be a whole number")
  }
  above <- matrix(10.3, nrow = 9, ncol = 4)
  expect_length(qc_chart(above, "xbar", center = 10, sigma = 1, run_length = Inf)$signals$run, 0)
  expect_equal(qc_chart(above, "xbar", center = 10, sigma = 1, run_length = 2)$signals$run, 2:9)
})

test_that("invalid center and sigma are errors", {
  for (bad in list(NaN, NA_real_, Inf, c(1, 2), "a")) {
    expect_error(qc_chart(x, "xbar", center = bad), "`center` must be")
    expect_error(qc_chart(x, "xbar", sigma = bad), "`sigma` must be")
  }
  expect_error(qc_chart(x, "xbar", sigma = 0), "`sigma` must be a single positive")
  expect_error(qc_chart(x, "xbar", sigma = -1), "`sigma` must be a single positive")
  # a negative center is fine for an X-bar chart but not for an R chart
  expect_no_error(qc_chart(x, "xbar", center = -5))
  expect_error(qc_chart(x, "R", center = -1), "`center` must be a single positive")
  expect_error(qc_chart(x, "R", center = 0), "`center` must be a single positive")
})

test_that("zero-variance subgroups give a warning", {
  flat <- matrix(5, 4, 3)
  expect_warning(q <- qc_chart(flat, "xbar"), "estimated sigma is 0")
  expect_equal(q$sigma, 0)
  expect_warning(qc_chart(flat, "R"), "estimated sigma is 0")
  # known sigma: no warning, even though the data have no spread
  expect_no_warning(qc_chart(flat, "xbar", sigma = 1))
})

test_that("a wide data frame error points to as.matrix()", {
  expect_error(qc_chart(as.data.frame(x), "xbar"), "as.matrix")
  expect_equal(qc_chart(as.matrix(as.data.frame(x)), "xbar")$center, qc_chart(x, "xbar")$center)
})

test_that("plot() draws the chart and returns the object invisibly", {
  q <- qc_chart(x, "xbar")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_invisible(res <- plot(q, title = "Plot method"))
  expect_identical(res, q)
})
