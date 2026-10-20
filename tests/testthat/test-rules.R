test_that("rule names are expanded and checked", {
  expect_equal(sqc:::expand_rules("multirule"), c("1-3s", "2-2s", "R-4s", "4-1s", "10x"))
  expect_equal(sqc:::expand_rules(c("1-2s", "multirule")), c("1-2s", "1-3s", "2-2s", "R-4s", "4-1s", "10x"))
  expect_error(sqc:::expand_rules("3-3s"), "Unknown rule")
  expect_error(apply_rules(c(0, 1), rules = character()), "character vector")
})

test_that("1-2s and 1-3s flag single controls", {
  z <- cbind(c(0, 2.4, 0, 3.2), c(0, 0, -2.1, 0))
  expect_equal(apply_rules(z, "1-2s")[, c("run", "level")], tibble::tibble(run = c(2L, 3L, 4L), level = c(1L, 2L, 1L)))
  v <- apply_rules(z, "1-3s")
  expect_equal(v$run, 4L)
  expect_equal(v$level, 1L)
})

test_that("2-2s needs two beyond the same side, within a run or across runs", {
  # across runs, one level
  expect_equal(apply_rules(c(0, 2.5, 2.6, 0), "2-2s")$run, 3L)
  # opposite sides do not count
  expect_equal(nrow(apply_rules(c(0, 2.5, -2.6, 0), "2-2s")), 0)
  # within a run, two levels beyond the same limit
  v <- apply_rules(cbind(c(0, 2.5, 0), c(0, 2.2, 0)), "2-2s")
  expect_equal(v$run, 2L)
  expect_true(is.na(v$level))
  # within a run, opposite sides are R-4s, not 2-2s
  expect_equal(nrow(apply_rules(cbind(c(2.5, 0), c(-2.5, 0)), "2-2s")), 0)
})

test_that("R-4s needs one control above +2 SD and another below -2 SD in the same run", {
  v <- apply_rules(cbind(c(0, 2.1, 0), c(0, -2.3, 0)), "R-4s")
  expect_equal(v$run, 2L)
  expect_true(is.na(v$level))
  expect_equal(nrow(apply_rules(cbind(c(2.1, 0), c(1.9, 0)), "R-4s")), 0)
  expect_equal(nrow(apply_rules(c(2.5, -2.5), "R-4s")), 0)   # one level: nothing to compare
})

test_that("4-1s and 10x need consecutive runs on one level", {
  expect_equal(apply_rules(c(0, 1.5, 1.2, 1.1, 1.3, 0), "4-1s")$run, 5L)
  expect_equal(nrow(apply_rules(c(1.5, 1.2, 0.9, 1.3, 1.1), "4-1s")), 0)
  expect_equal(apply_rules(c(-1.5, -1.2, -1.1, -1.3), "4-1s")$run, 4L)
  expect_equal(apply_rules(rep(0.3, 10), "10x")$run, 10L)
  expect_equal(nrow(apply_rules(rep(0.3, 9), "10x")), 0)
  expect_equal(nrow(apply_rules(c(rep(0.3, 5), -0.1, rep(0.3, 4)), "10x")), 0)
  # the two levels are counted separately: alternating levels do not make a run
  z <- cbind(rep(c(1.5, 0), 5), rep(c(0, 1.5), 5))
  expect_equal(nrow(apply_rules(z, "4-1s")), 0)
})

test_that("apply_rules returns an empty tibble when nothing fires, and validates z", {
  v <- apply_rules(c(0.1, -0.2, 0.3))
  expect_s3_class(v, "tbl_df")
  expect_equal(names(v), c("run", "rule", "level"))
  expect_equal(nrow(v), 0)
  expect_error(apply_rules(c(0, NA)), "no missing")
})

test_that("null rejection probabilities per point match the derivations", {
  set.seed(11)
  n <- 2e6
  z1 <- matrix(rnorm(n), ncol = 1)                      # one level
  f1 <- sqc:::rule_fires(z1, c("1-3s", "2-2s", "4-1s", "10x"))
  q <- function(k) 1 - pnorm(k)
  expected <- c("1-3s" = 2 * q(3), "2-2s" = 2 * q(2)^2, "4-1s" = 2 * q(1)^4, "10x" = 2 * 0.5^10)
  expect_equal(colMeans(f1), expected, tolerance = 0.15)   # within 15 percent (Monte Carlo)
  # R-4s with two levels
  z2 <- matrix(rnorm(n), ncol = 2)
  r4 <- mean(sqc:::rule_fires(z2, "R-4s"))
  expect_equal(r4, 2 * q(2)^2, tolerance = 0.1)
  # the weaker reading "range between the two controls exceeds 4 SD" has a different probability
  expect_equal(mean(abs(z2[, 1] - z2[, 2]) > 4), 2 * (1 - pnorm(4 / sqrt(2))), tolerance = 0.1)
})
