test_that("T-squared with known parameters is the standardized distance and has the chi-square limit", {
  S <- matrix(c(36, 0.8 * 6 * 4, 0.8 * 6 * 4, 16), 2)
  d <- data.frame(a = c(125, 138, 118), b = c(80, 70, 83))
  h <- hotelling_chart(d, mean = c(125, 80), cov = S)
  z1 <- (d$a - 125) / 6; z2 <- (d$b - 80) / 4
  expect_equal(h$statistics$t2, (z1^2 - 2 * 0.8 * z1 * z2 + z2^2) / (1 - 0.8^2))
  expect_equal(h$statistics$t2[1], 0)
  expect_equal(h$ucl[["known"]], qchisq(1 - 0.0027, 2)); expect_true(h$known); expect_true(is.na(h$m))
  expect_equal(h$signals, 2L)                            # 138/70 is far outside the ellipse
  expect_equal(unique(h$statistics$phase), "known")
  # the limit follows alpha
  expect_equal(hotelling_chart(d, mean = c(125, 80), cov = S, alpha = 0.01)$ucl[["known"]], qchisq(0.99, 2))
})

test_that("estimated parameters use the beta limit in Phase I and the F limit in Phase II", {
  h <- hotelling_chart(bp_home, c("systolic", "diastolic"), phase1 = 1:60)
  m <- 60; p <- 2
  expect_equal(h$ucl[["I"]], (m - 1)^2 / m * qbeta(1 - 0.0027, 1, (m - 3) / 2))
  expect_equal(h$ucl[["II"]], p * (m + 1) * (m - 1) / (m * (m - p)) * qf(1 - 0.0027, p, m - p))
  expect_gt(h$ucl[["II"]], qchisq(1 - 0.0027, 2))        # estimation widens the limit
  expect_equal(h$statistics$phase, rep(c("I", "II"), c(60, 30)))
  expect_equal(h$statistics$ucl, rep(unname(h$ucl), c(60, 30)))
  expect_equal(h$mean, colMeans(bp_home[1:60, c("systolic", "diastolic")]))
  expect_equal(unname(h$cov), unname(cov(bp_home[1:60, c("systolic", "diastolic")])))
  expect_equal(h$statistics$t2, mahalanobis(as.matrix(bp_home[c("systolic", "diastolic")]), h$mean, h$cov))
  expect_equal(h$m, 60)
  expect_output(print(h), "Phase II")
  expect_s3_class(as.data.frame(h), "data.frame")
  # a logical phase1 gives the same chart
  h2 <- hotelling_chart(bp_home, c("systolic", "diastolic"), phase1 = bp_home$day <= 60)
  expect_equal(h2$statistics$t2, h$statistics$t2)
})

test_that("the blood pressure example: 138/70 is inside the box but outside the ellipse, 144/91 the reverse", {
  h <- hotelling_chart(bp_home, c("systolic", "diastolic"), phase1 = 1:60)
  expect_equal(h$signals, 75L)
  for (v in c("systolic", "diastolic")) {
    x <- bp_home[[v]]; mr <- mean(abs(diff(x[1:60])))
    lo <- mean(x[1:60]) - 2.66 * mr; hi <- mean(x[1:60]) + 2.66 * mr
    if (v == "systolic") { expect_true(x[75] > lo && x[75] < hi); expect_true(x[82] > hi); expect_equal(which(x < lo | x > hi), 82L) } else
      expect_equal(which(x < lo | x > hi), integer(0))
  }
  expect_lt(h$statistics$t2[82], h$ucl[["II"]])
  expect_gt(h$statistics$t2[75], 2 * h$ucl[["II"]])
})

test_that("contributions measure the extra distance of each variable given the others", {
  h <- hotelling_chart(bp_home, c("systolic", "diastolic"), phase1 = 1:60)
  t2 <- h$statistics$t2; z <- as.matrix(h$z); cn <- as.matrix(h$contributions)
  expect_equal(cn[, "systolic"], t2 - z[, "diastolic"]^2, ignore_attr = TRUE)       # p = 2: T2 minus the other variable's squared z
  expect_equal(cn[, "diastolic"], t2 - z[, "systolic"]^2, ignore_attr = TRUE)
  expect_gt(cn[75, "diastolic"], cn[75, "systolic"])                              # the diastolic is unusual given the systolic
  expect_true(all(cn >= -1e-8))
  # three variables: the contribution is T2 minus the T2 of the other two
  set.seed(2)
  d3 <- data.frame(a = rnorm(40), b = rnorm(40), c = rnorm(40)); d3$b <- d3$b + d3$a
  h3 <- hotelling_chart(d3)
  expect_equal(h3$p, 3)
  idx <- 7; other <- c("a", "c")
  t2_ac <- mahalanobis(as.matrix(d3[other]), colMeans(d3[other]), cov(d3[other]))
  expect_equal(h3$contributions$b[idx], h3$statistics$t2[idx] - t2_ac[idx])
  expect_equal(h3$vars, c("a", "b", "c"))
})

test_that("t2_ellipse traces the points of constant T-squared", {
  S <- matrix(c(4, 3, 3, 9), 2); mu <- c(10, -2)
  e <- t2_ellipse(mu, S, c = 5.99)
  expect_equal(mahalanobis(as.matrix(e), mu, S), rep(5.99, nrow(e)))
  expect_equal(nrow(t2_ellipse(mu, S, n = 50)), 50)
  expect_error(t2_ellipse(c(1, 2, 3), S)); expect_error(t2_ellipse(mu, matrix(1, 2, 2))); expect_error(t2_ellipse(mu, S, c = -1))
  p <- autoplot(hotelling_chart(bp_home, c("systolic", "diastolic"), phase1 = 1:60), which = "ellipse")
  expect_s3_class(p, "ggplot")
  expect_s3_class(autoplot(hotelling_chart(bp_home, c("systolic", "diastolic"), phase1 = 1:60)), "ggplot")
  expect_error(autoplot(hotelling_chart(data.frame(a = rnorm(20), b = rnorm(20), c = rnorm(20))), which = "ellipse"))
})

test_that("the box has about twice the false alarm probability of the ellipse, and misses what the ellipse sees", {
  set.seed(4)
  s0 <- sim_box_ellipse(0, n = 300000)
  expect_equal(s0$box, 1 - (1 - 0.0027)^2, tolerance = 0.06); expect_equal(s0$ellipse, 0.0027, tolerance = 0.1)
  s8 <- sim_box_ellipse(0.8, n = 300000)
  expect_equal(s8$ellipse, 0.0027, tolerance = 0.1)
  expect_lt(s8$box, s0$box)                                      # correlated charts alarm together
  expect_equal(s8$box_only + s8$ellipse_only + s8$both + s8$neither, 1)
  expect_gt(s8$ellipse_only / s8$ellipse, s0$ellipse_only / s0$ellipse)   # more of the ellipse's signals are inside the box as rho grows
  expect_gt(s8$box_only, s8$ellipse_only)                                   # in control, the box has the extra false alarms
  # the variables stop moving together: the box sees almost nothing, the ellipse sees it
  sd <- sim_box_ellipse(0.8, rho_true = 0)
  expect_lt(sd$box, 0.02); expect_gt(sd$ellipse, 0.1)
  # a shift of one variable
  ss <- sim_box_ellipse(0.8, shift = c(0, 1.5))
  expect_gt(ss$ellipse, 2 * ss$box)
  s00 <- sim_box_ellipse(0, shift = c(0, 1.5))                     # uncorrelated variables: the ellipse has no advantage
  expect_gt(s00$box, s00$ellipse)
  p <- sim_box_ellipse(0.8, n = 500, points = TRUE)
  expect_equal(nrow(p), 500); expect_setequal(names(p), c("x", "y", "t2", "box", "ellipse", "category"))
  expect_error(sim_box_ellipse(1)); expect_error(sim_box_ellipse(0.5, shift = 1)); expect_error(sim_box_ellipse(0.5, n = 3))
})

test_that("with estimated parameters the chi-square limit is wrong for small m and the F and beta limits are right", {
  set.seed(6)
  s <- sim_hotelling(m = 20, reps = 1000)
  r <- function(ph, l) s$rate[s$phase == ph & s$limit == l]
  expect_equal(r("II", "F"), 0.0027, tolerance = 0.3)
  expect_equal(r("I", "beta"), 0.0027, tolerance = 0.4)
  expect_gt(r("II", "chisq"), 2.5 * r("II", "F"))
  expect_lt(r("I", "chisq"), r("I", "beta"))
  expect_equal(s$ucl[s$limit == "F"], 2 * 21 * 19 / (20 * 18) * qf(1 - 0.0027, 2, 18))
  expect_equal(nrow(s), 4)
  expect_error(sim_hotelling(m = 3)); expect_error(sim_hotelling(m = 20, alpha = 2))
})

test_that("hotelling_chart checks its arguments", {
  d <- data.frame(a = rnorm(30), b = rnorm(30), g = letters[rep(1:3, 10)])
  expect_error(hotelling_chart(d, "a")); expect_error(hotelling_chart(d, c("a", "zz")))
  expect_s3_class(hotelling_chart(d), "hotelling_chart")                                # the numeric columns are used
  expect_equal(hotelling_chart(d)$vars, c("a", "b"))
  expect_error(hotelling_chart(data.frame(a = c(1:29, NA), b = rnorm(30))))
  expect_error(hotelling_chart(d, c("a", "b"), phase1 = 1:3))                             # too few rows
  expect_error(hotelling_chart(d, c("a", "b"), phase1 = 0:10))
  expect_error(hotelling_chart(d, c("a", "b"), mean = c(0, 0)))                           # both or neither
  expect_error(hotelling_chart(d, c("a", "b"), mean = c(0, 0), cov = diag(3)))
  expect_error(hotelling_chart(d, c("a", "b"), mean = c(0, 0), cov = matrix(c(1, 1, 1, 1), 2)))
  e <- data.frame(a = rnorm(30)); e$b <- 2 * e$a
  expect_error(hotelling_chart(e))                                                         # singular covariance
  expect_error(hotelling_chart(d, c("a", "b"), alpha = 0))
})
