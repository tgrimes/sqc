test_that("1-3s and 1-2s ARLs match the exact geometric values", {
  set.seed(1)
  a13 <- sim_run_length("1-3s", reps = 2000)
  p13 <- 1 - (1 - 2 * (1 - pnorm(3)))^2
  se13 <- sd(a13) / sqrt(2000)
  expect_lt(abs(arl(a13) - 1 / p13), 5 * se13)
  a12 <- sim_run_length("1-2s", reps = 3000)
  p12 <- 1 - (1 - 2 * (1 - pnorm(2)))^2
  expect_lt(abs(arl(a12) - 1 / p12), 5 * sd(a12) / sqrt(3000))
})

test_that("pattern rules do not have ARL = 1/p: 2-2s and 10x match the exact values", {
  set.seed(2)
  q <- 1 - pnorm(2)
  exact_22 <- (1 + q) / (2 * q^2)           # Markov chain result for one level
  a22 <- sim_run_length("2-2s", levels = 1, reps = 2500)
  expect_lt(abs(arl(a22) - exact_22), 5 * sd(a22) / sqrt(2500))
  expect_gt(exact_22, 1 / (2 * q^2))         # longer than 1/p
  a10 <- sim_run_length("10x", levels = 1, reps = 1500)
  expect_lt(abs(arl(a10) - (2^10 - 1)), 5 * sd(a10) / sqrt(1500))   # 1023, twice 1/p = 512
})

test_that("shifts shorten the run length and the multirule beats 1-3s", {
  set.seed(3)
  in_control <- arl(sim_run_length("multirule", reps = 600))
  shifted <- arl(sim_run_length("multirule", shift = 2, reps = 600))
  expect_lt(shifted, 4)
  expect_lt(shifted, in_control)
  expect_lt(arl(sim_run_length("multirule", shift = 1, reps = 600)), arl(sim_run_length("1-3s", shift = 1, reps = 600)))
})

test_that("lognormal data have a much shorter in-control ARL for the 1-3s rule", {
  set.seed(4)
  rl <- sim_run_length("1-3s", levels = 1, dist = "lognormal", sdlog = 0.5, reps = 2500)
  s <- 0.5
  m <- exp(s^2 / 2); sd_ <- sqrt((exp(s^2) - 1) * exp(s^2))
  p <- 1 - plnorm(m + 3 * sd_, 0, s)                    # the lower limit is below 0
  expect_lt(abs(arl(rl) - 1 / p), 5 * sd(rl) / sqrt(2500))
  expect_lt(arl(rl), 100)
})

test_that("results print, can be summarized, and report censoring", {
  set.seed(5)
  rl <- sim_run_length("1-3s", reps = 20)
  expect_s3_class(rl, "sqc_run_lengths")
  expect_output(print(rl), "ARL")
  expect_length(rl, 20)
  cens <- sim_run_length("10x", levels = 1, reps = 5, max_runs = 20)
  expect_output(print(cens), "censored")
  expect_equal(max(unclass(cens)), 20)
})

test_that("arl_table returns a tidy table", {
  set.seed(6)
  tab <- arl_table(list("1-3s", Multirule = "multirule"), shifts = c(0, 2), reps = 200)
  expect_equal(names(tab), c("rule", "shift", "arl", "sd", "se"))
  expect_equal(tab$rule, c("1-3s", "1-3s", "Multirule", "Multirule"))
  expect_equal(tab$shift, c(0, 2, 0, 2))
  expect_true(all(tab$arl > 0))
  tab2 <- arl_table(c("1-3s", "1-2s"), shifts = 0, reps = 100)
  expect_equal(tab2$rule, c("1-3s", "1-2s"))
})

test_that("arguments are validated", {
  expect_error(sim_run_length("1-3s", levels = 0), "`levels` must be")
  expect_error(sim_run_length("1-3s", reps = 2.5), "`reps` must be")
  expect_error(sim_run_length("1-3s", shift = c(1, 2, 3)), "`shift` must be")
  expect_error(sim_run_length("1-3s", sdlog = -1), "`sdlog` must be")
  expect_error(sim_run_length("R-4s", levels = 1), "R-4s rule needs")
  expect_error(arl(1:5), "must come from")
  expect_error(sim_run_length("nonsense"), "Unknown rule")
})


test_that("scale changes the spread of the standardized results", {
  set.seed(21)
  # one level, 1-3s: a point is beyond 3 with probability 2 * (1 - pnorm(3 / scale))
  for (sc in c(0.9, 1.5)) {
    p <- 2 * (1 - pnorm(3 / sc))
    rl <- sim_run_length("1-3s", levels = 1, scale = sc, reps = 4000)
    expect_lt(abs(arl(rl) - 1 / p), 4.5 * stats::sd(unclass(rl)) / sqrt(4000))
  }
  # shift and scale together describe estimated limits: z = (X - mu_hat) / sigma_hat
  set.seed(22)
  mu_hat <- 0.4; sigma_hat <- 0.8                     # true mean 0 and SD 1
  p <- 1 - pnorm(mu_hat + 3 * sigma_hat) + pnorm(mu_hat - 3 * sigma_hat)
  rl <- sim_run_length("1-3s", levels = 1, shift = -mu_hat / sigma_hat, scale = 1 / sigma_hat, reps = 4000)
  expect_lt(abs(arl(rl) - 1 / p), 4.5 * stats::sd(unclass(rl)) / sqrt(4000))
  expect_equal(attr(rl, "scale"), 1 / sigma_hat)
})

test_that("scale is validated and one value per level is allowed", {
  expect_error(sim_run_length("1-3s", scale = 0), "`scale` must be a positive")
  expect_error(sim_run_length("1-3s", scale = c(1, 2, 3)), "`scale` must be a positive")
  expect_error(sim_run_length("1-3s", scale = NA_real_), "`scale` must be a positive")
  set.seed(23)
  rl <- sim_run_length("1-3s", levels = 2, scale = c(1, 1), reps = 50)
  expect_s3_class(rl, "sqc_run_lengths")
})
