set.seed(21)
x <- cbind(rnorm(30, 100, 3), rnorm(30, 250, 6))
x[19:30, ] <- x[19:30, ] + rep(c(4.5, 9), each = 12)     # a 1.5 SD shift from run 19

test_that("lj_chart standardizes and applies the rules", {
  q <- lj_chart(x, mean = c(100, 250), sd = c(3, 6))
  expect_s3_class(q, "lj_chart")
  expect_equal(q$z[, 1], (x[, 1] - 100) / 3)
  expect_equal(q$z[, 2], (x[, 2] - 250) / 6)
  expect_equal(q$violations, apply_rules(q$z))
  expect_equal(q$rejected, sort(unique(q$violations$run)))
  expect_true(all(q$rejected >= 19))          # nothing rejected before the shift with this seed
})

test_that("scalar mean and sd are recycled and a vector means one level", {
  q <- lj_chart(rnorm(10), mean = 0, sd = 1)
  expect_equal(q$n_levels, 1)
  q2 <- lj_chart(cbind(rnorm(5), rnorm(5)), mean = 0, sd = 1, rules = "1-2s")
  expect_equal(q2$mean, c(0, 0))
  expect_equal(q2$rules, "1-2s")
})

test_that("print, as.data.frame, and autoplot work", {
  q <- lj_chart(x, mean = c(100, 250), sd = c(3, 6))
  expect_output(print(q), "30 runs")
  expect_output(print(q), "runs rejected")
  d <- as.data.frame(q)
  expect_equal(nrow(d), 60)
  expect_equal(names(d), c("run", "level", "value", "z", "rejected"))
  p <- ggplot2::autoplot(q, title = "Test", subtitle = "sub")
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  expect_equal(p$labels$title, "Test")
  expect_no_error(ggplot2::ggplot_build(ggplot2::autoplot(lj_chart(rnorm(10), 0, 1))))   # nothing rejected
})

test_that("labels can be set and bad input gives clear errors", {
  q <- lj_chart(x, mean = c(100, 250), sd = c(3, 6), level_labels = c("Low", "High"))
  expect_equal(q$level_labels, c("Low", "High"))
  expect_error(lj_chart(x, mean = c(100, 250), sd = c(3, 6), level_labels = "Low"), "one label per level")
  expect_error(lj_chart(x, mean = 100:102, sd = 3), "`mean` must be")
  expect_error(lj_chart(x, mean = 100, sd = -1), "`sd` must be positive")
  expect_error(lj_chart(cbind(c(1, NA), c(1, 2)), 0, 1), "no missing")
  expect_error(lj_chart(x, 100, 3, rules = "foo"), "Unknown rule")
})
