test_that("sim_ar1 has the requested marginal standard deviation and autocorrelations", {
  set.seed(1)
  x <- sim_ar1(40000, phi = 0.8, sigma = 2, mean = 10)
  expect_equal(mean(x), 10, tolerance = 0.08)
  expect_equal(sd(x), 2, tolerance = 0.03)
  a <- acf(x, 3, plot = FALSE)$acf[2:4]
  expect_equal(a, 0.8^(1:3), tolerance = 0.03)
  y <- sim_ar1(1000, phi = 0, sigma = 1, shift = 5, change = 501)
  expect_gt(mean(y[501:1000]) - mean(y[1:500]), 4.5)
  expect_equal(length(y), 1000)
  expect_error(sim_ar1(100, phi = 1)); expect_error(sim_ar1(100, phi = 0.5, sigma = -1)); expect_error(sim_ar1(100, phi = 0.5, shift = 1))
})

test_that("the false alarm probability of the I chart follows the moving range bias", {
  expect_equal(i_chart_alarm_rate(0), 2 * pnorm(-3))
  expect_equal(i_chart_alarm_rate(0.5), 2 * pnorm(-3 * sqrt(0.5)))
  expect_equal(i_chart_alarm_rate(0.9), 0.343, tolerance = 0.003)
  expect_equal(i_chart_alarm_rate(c(0, 0.5, 0.9), k = 2), 2 * pnorm(-2 * sqrt(1 - c(0, 0.5, 0.9))))
  expect_error(i_chart_alarm_rate(1)); expect_error(i_chart_alarm_rate(0.5, k = -1))
  # the moving range estimates sigma * sqrt(1 - rho)
  set.seed(3)
  x <- sim_ar1(50000, phi = 0.7)
  expect_equal(mean(abs(diff(x))) / 1.128, sqrt(1 - 0.7), tolerance = 0.03)
  expect_equal(mean(diff(x)^2), 2 * (1 - 0.7), tolerance = 0.03)
})

test_that("the simulated false alarm rates agree with the theory", {
  set.seed(5)
  s0 <- sim_dependent(0, reps = 300); s5 <- sim_dependent(0.5, reps = 300); s9 <- sim_dependent(0.9, reps = 300)
  m <- function(s, ch) mean(s$share[s$chart == ch])
  expect_equal(m(s0, "I"), 0.0027, tolerance = 0.5)
  expect_equal(m(s5, "I"), i_chart_alarm_rate(0.5), tolerance = 0.25)
  expect_equal(m(s9, "I"), i_chart_alarm_rate(0.9), tolerance = 0.2)
  expect_lt(m(s9, "residual"), 0.006)                                   # the residual chart keeps the nominal rate
  expect_lt(m(s5, "residual"), 0.006)
  expect_gt(m(s9, "I"), 50 * m(s9, "residual"))
  expect_equal(nrow(s9), 900); expect_setequal(unique(s9$chart), c("I", "I_sd", "residual"))
  expect_error(sim_dependent(1)); expect_error(sim_dependent(0.5, n = 3))
})

test_that("the residual chart uses the one-step-ahead residuals", {
  x <- c(5, 5.4, 5.1, 4.8, 5.2, 5.0, 5.3, 4.9, 5.1, 5.6, 5.0)
  a <- ar1_chart(x, phi = 0.6, mean = 5, sigma_e = 0.2, k = 3)
  expect_true(is.na(a$statistics$residual[1]))
  expect_equal(a$statistics$fitted[-1], 5 + 0.6 * (x[-11] - 5))
  expect_equal(a$statistics$residual[-1], x[-1] - (5 + 0.6 * (x[-11] - 5)))
  expect_equal(a$statistics$ucl, rep(0.6, 11)); expect_equal(a$sigma_x, 0.2 / sqrt(1 - 0.36))
  expect_equal(a$signals, which(abs(a$statistics$residual) > 0.6))
  expect_output(print(a), "AR\\(1\\) fit")
  expect_s3_class(autoplot(a), "ggplot"); expect_s3_class(autoplot(a, which = "data"), "ggplot"); expect_s3_class(as.data.frame(a), "data.frame")
  expect_error(ar1_chart(x[1:5])); expect_error(ar1_chart(x, phi = 1.2, mean = 5)); expect_error(ar1_chart(c(x, NA)))
})

test_that("the fitted parameters are close to the truth and the residuals are independent", {
  set.seed(9)
  x <- sim_ar1(2000, phi = 0.85, sigma = 1.5, mean = 20)
  a <- ar1_chart(x)
  expect_equal(a$phi, 0.85, tolerance = 0.05); expect_equal(a$mean, 20, tolerance = 0.3)
  expect_equal(a$sigma_e, 1.5 * sqrt(1 - 0.85^2), tolerance = 0.06)
  expect_lt(abs(cor(a$statistics$residual[-1][-1], a$statistics$residual[-1][-1999])), 0.08)
  expect_lt(length(a$signals), 12)
  b <- ar1_chart(x, sigma_from = "sd"); expect_equal(b$sigma_e, sqrt(mean(a$statistics$residual[-1]^2)))
})

test_that("the refrigerator data show the story", {
  x <- fridge_temp$temp_c
  expect_equal(nrow(fridge_temp), 960); expect_equal(as.numeric(diff(fridge_temp$time[1:2]), units = "mins"), 15)
  expect_true(all(x > 2 & x < 8))
  expect_equal(cor(x[-1], x[-960]), 0.9, tolerance = 0.02)
  q <- qc_chart(fridge_temp, "I", value = temp_c)
  expect_gt(length(q$signals$beyond) / 960, 0.25)                          # the I chart alarms for a third of the readings
  a <- ar1_chart(x, k = 4)
  expect_equal(a$signals, 729L)
  expect_equal(format(fridge_temp$time[729], "%Y-%m-%d %H:%M"), "2026-01-12 14:00")
  expect_gt(length(ar1_chart(x, k = 3)$signals), 1)
  sub <- x[seq(1, 960, by = 20)]
  expect_equal(length(qc_chart(data.frame(v = sub), "I", value = v)$signals$beyond), 0)
  expect_s3_class(acf_plot(x, lag_max = 40), "ggplot")
})

test_that("acf_plot returns the autocorrelations", {
  set.seed(2)
  x <- sim_ar1(300, 0.7)
  p <- acf_plot(x, lag_max = 10)
  expect_equal(p$data$lag, 1:10)
  expect_equal(p$data$acf, as.numeric(acf(x, 10, plot = FALSE)$acf)[-1])
  expect_error(acf_plot(1:3)); expect_error(acf_plot(c(1:10, NA)))
})
