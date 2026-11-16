test_that("the Poisson CUSUM uses the observed-minus-expected log-likelihood ratio", {
  obs <- c(3, 5, 2, 6, 7, 4)
  ex  <- c(4, 4.5, 3, 4, 4.2, 4)
  cs <- cusum_poisson(obs, ex, rr = 1.5, h = 2)
  w <- obs * log(1.5) - ex * 0.5
  expect_equal(cs$statistics$weight, w)
  s <- numeric(length(w)); prev <- 0
  for (t in seq_along(w)) { prev <- max(0, prev + w[t]); s[t] <- prev }
  expect_equal(cs$statistics$statistic, s)
  expect_equal(cs$statistics$signal, s > 2)
  # a Poisson count equal to its expectation has a negative weight for a relative risk above 1
  expect_lt(cusum_poisson(4, 4, rr = 1.5)$statistics$weight, 0)
  # the weight is the log-likelihood ratio: positive once the count is large enough relative to the expected count
  expect_equal(cusum_poisson(6, 4, rr = 1.5)$statistics$weight, 6 * log(1.5) - 4 * 0.5)
  expect_gt(cusum_poisson(6, 4, rr = 1.5)$statistics$weight, 0)
})

test_that("the signal reports the estimated start and the accumulated excess", {
  obs <- c(4, 5, 4, 3, 7, 8, 9, 8, 10)
  ex  <- rep(4, 9)
  cs <- cusum_poisson(obs, ex, rr = 1.5, h = 2)
  s <- cs$signals[1, ]
  st <- cs$statistics$statistic
  expect_equal(s$index, cs$first_signal)
  expect_equal(s$start, max(which(st[seq_len(s$index)] == 0)) + 1)
  expect_equal(s$excess, sum(obs[s$start:s$index] - ex[s$start:s$index]))
  expect_output(print(cs), "estimated to start")
  expect_s3_class(as.data.frame(cs), "data.frame")
  expect_no_error(autoplot(cs))
})

test_that("a reset restarts the sum after each signal", {
  obs <- rep(10, 12); ex <- rep(4, 12)
  once <- cusum_poisson(obs, ex, rr = 1.5, h = 5)
  more <- cusum_poisson(obs, ex, rr = 1.5, h = 5, reset = TRUE)
  expect_equal(once$first_signal, more$first_signal)
  expect_gt(nrow(more$signals), 1)
  expect_lt(more$statistics$statistic[more$first_signal + 1], once$statistics$statistic[once$first_signal + 1])
})

test_that("the risk-adjusted CUSUM uses the Bernoulli log-likelihood ratio", {
  y <- c(0, 1, 0, 0, 1, 1)
  p <- c(0.05, 0.1, 0.2, 0.3, 0.1, 0.4)
  cs <- cusum_risk(y, p, odds_ratio = 2, h = 1)
  w <- y * log(2) - log(1 - p + 2 * p)
  expect_equal(cs$statistics$weight, w)
  # a death adds more weight when the predicted risk is low than when it is high
  expect_gt(cusum_risk(1, 0.01, 2)$statistics$weight, cusum_risk(1, 0.5, 2)$statistics$weight)
  # the weights have mean zero under the null of odds ratio 1 only for the matching alternative: the mean weight is negative
  p0 <- 0.1; n0 <- 1 - p0 + 2 * p0
  expect_lt(p0 * (log(2) - log(n0)) + (1 - p0) * (-log(n0)), 0)
  # an odds ratio below 1 detects improvement: weights are positive for survivals' complement (deaths get negative weight)
  expect_lt(cusum_risk(1, 0.1, 0.5)$statistics$weight, 0)
})

test_that("the expected number of events is the sum of the risks", {
  d <- surgeon_cases[surgeon_cases$surgeon == "C", ]
  cs <- cusum_risk(d$death, d$risk, 2, h = 4)
  expect_equal(cs$statistics$expected, d$risk)
  if (!is.na(cs$first_signal)) {
    s <- cs$signals[1, ]
    expect_equal(s$excess, sum(d$death[s$start:s$index] - d$risk[s$start:s$index]))
  }
})

test_that("the CUSUMs check their arguments", {
  expect_error(cusum_poisson(c(1, 2), c(1, 2, 3)))
  expect_error(cusum_poisson(c(1.5, 2), c(1, 2)))
  expect_error(cusum_poisson(c(1, -2), c(1, 2)))
  expect_error(cusum_poisson(c(1, 2), c(0, 2)))
  expect_error(cusum_poisson(c(1, 2), c(1, 2), rr = 1))
  expect_error(cusum_risk(c(0, 2), c(0.1, 0.2)))
  expect_error(cusum_risk(c(0, 1), c(0.1, 1)))
  expect_error(cusum_risk(c(0, 1), c(0, 0.5)))
  expect_error(cusum_risk(c(0, 1), c(0.1, 0.2), odds_ratio = 1))
  expect_error(cusum_risk(c(0, 1), c(0.1, 0.2), h = -1))
})

test_that("the simulated maximum gives a false alarm probability", {
  set.seed(5)
  ex <- gp_deaths$expected
  s <- sim_cusum_poisson(ex, rr = 1.5, reps = 4000)
  expect_equal(nrow(s), 4000)
  h <- cusum_limit(s, alpha = 0.05)
  expect_equal(mean(s$max_stat > h), 0.05, tolerance = 0.2)
  expect_lte(mean(s$max_stat > h), 0.05)
  expect_gt(cusum_limit(s, alpha = 0.01), h)
  # with the matching shift the chart signals much sooner than in control
  set.seed(6)
  late <- sim_cusum_poisson(ex, rr = 1.5, true_rr = 1.5, h = h, reps = 2000)
  expect_gt(mean(!is.na(late$run_length)), 0.8)
  expect_lt(median(late$run_length, na.rm = TRUE), length(ex))
  expect_error(cusum_limit(data.frame(a = 1)))
  expect_error(cusum_limit(s, alpha = 2))
})

test_that("the risk-adjusted simulation is calibrated on the predicted risks", {
  set.seed(8)
  r <- surgeon_cases$risk[surgeon_cases$surgeon == "A"]
  s <- sim_cusum_risk(r, odds_ratio = 2, reps = 3000)
  h <- cusum_limit(s, 0.05)
  expect_equal(mean(s$max_stat > h), 0.05, tolerance = 0.25)
  bad <- sim_cusum_risk(r, odds_ratio = 2, true_odds_ratio = 2, h = h, reps = 1500)
  expect_gt(mean(!is.na(bad$run_length)), mean(s$max_stat > h))
  expect_error(sim_cusum_risk(c(0.1, 1.2)))
  expect_error(sim_cusum_risk(c(0.1, 0.2), odds_ratio = 1))
  expect_error(sim_cusum_poisson(c(1, 0)))
})

test_that("the simulated datasets tell the stories they are built for", {
  expect_equal(dim(gp_deaths), c(24, 3))
  expect_true(all(abs((gp_deaths$observed - gp_deaths$expected) / sqrt(gp_deaths$expected)) < 3))     # no Shewhart signal
  p <- cusum_poisson(gp_deaths$observed, gp_deaths$expected, rr = 1.5, h = 4.28)
  expect_false(is.na(p$first_signal))
  expect_lt(p$first_signal, 15)
  expect_equal(dim(surgeon_cases), c(1200, 4))
  expect_true(all(surgeon_cases$risk > 0 & surgeon_cases$risk < 1))
  by_surgeon <- split(surgeon_cases, surgeon_cases$surgeon)
  adjusted <- vapply(by_surgeon, function(d) max(cusum_risk(d$death, d$risk, 2)$statistics$statistic), numeric(1))
  expect_true(adjusted[["C"]] > 4.5)                         # past the limit for a 5% false alarm probability (about 4.2)
  expect_true(all(adjusted[c("A", "B")] < 3.5))              # the other two surgeons stay below it
})
