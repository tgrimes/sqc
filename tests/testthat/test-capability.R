test_that("capability reproduces the textbook formulas", {
  x <- tablet_supplier$weight_mg
  cap <- capability(x, lsl = 475, usl = 525)
  m <- mean(x); s <- sd(x); n <- length(x)
  i <- setNames(cap$indices$estimate, cap$indices$index)
  expect_equal(i[["Cp"]], 50 / (6 * s))
  expect_equal(i[["Cpl"]], (m - 475) / (3 * s))
  expect_equal(i[["Cpu"]], (525 - m) / (3 * s))
  expect_equal(i[["Cpk"]], min((m - 475), (525 - m)) / (3 * s))
  expect_equal(i[["Cpk"]], 1.40, tolerance = 0.005)
  expect_equal(i[["Ppk"]], i[["Cpk"]])                          # no separate within sigma: the same
  expect_false("Cpm" %in% cap$indices$index)
  expect_equal(cap$n, 30)
})

test_that("the confidence intervals have the documented form", {
  x <- tablet_supplier$weight_mg
  cap <- capability(x, lsl = 475, usl = 525, conf = 0.9)
  i <- setNames(cap$indices$estimate, cap$indices$index)
  lo <- setNames(cap$indices$lower, cap$indices$index); up <- setNames(cap$indices$upper, cap$indices$index)
  n <- 30; z <- qnorm(0.95)
  se <- sqrt(1 / (9 * n) + i[["Cpk"]]^2 / (2 * (n - 1)))
  expect_equal(lo[["Cpk"]], i[["Cpk"]] - z * se); expect_equal(up[["Cpk"]], i[["Cpk"]] + z * se)
  expect_equal(lo[["Cp"]], i[["Cp"]] * sqrt(qchisq(0.05, n - 1) / (n - 1)))
  expect_equal(up[["Cp"]], i[["Cp"]] * sqrt(qchisq(0.95, n - 1) / (n - 1)))
  expect_equal(cap$lower_one_sided, i[["Cpk"]] - qnorm(0.9) * se)
  expect_lt(lo[["Cpk"]], i[["Cpk"]]); expect_gt(up[["Cpk"]], i[["Cpk"]])
  none <- capability(x, lsl = 475, usl = 525, ci = "none")
  expect_true(all(is.na(none$indices$lower)))
})

test_that("the observed and expected fractions outside the specification are computed", {
  x <- tablet_press$weight_mg
  cap <- capability(x, lsl = 475, usl = 525)
  expect_equal(cap$observed_out, 0.02)
  expect_equal(cap$expected_out, pnorm(475, mean(x), sd(x)) + pnorm(525, mean(x), sd(x), lower.tail = FALSE))
  expect_equal(cap$dpmo, 1e6 * cap$expected_out)
  expect_output(print(cap), "Outside the specification")
})

test_that("a supplied within-subgroup sigma separates the Cp indices from the Pp indices", {
  x <- tablet_press$weight_mg
  q <- qc_chart(tablet_press, "xbar", value = weight_mg, sample = subgroup)
  cap <- capability(x, lsl = 475, usl = 525, sigma = q$sigma, df = 20 * 4)
  i <- setNames(cap$indices$estimate, cap$indices$index)
  expect_equal(i[["Cp"]], 50 / (6 * q$sigma))
  expect_equal(i[["Pp"]], 50 / (6 * sd(x)))
  expect_false(is.na(cap$indices$lower[cap$indices$index == "Cp"]))
  # without df the intervals are not available
  expect_true(all(is.na(capability(x, lsl = 475, usl = 525, sigma = q$sigma)$indices$lower)))
})

test_that("one-sided specifications and the target index work", {
  x <- tablet_supplier$weight_mg
  up <- capability(x, usl = 525)
  expect_equal(up$indices$index, c("Cpu", "Cpk", "Ppk"))
  expect_equal(up$indices$estimate[1], (525 - mean(x)) / (3 * sd(x)))
  cp <- capability(x, lsl = 475, usl = 525, target = 500)
  expect_equal(cp$indices$estimate[cp$indices$index == "Cpm"], 50 / (6 * sqrt(var(x) + (mean(x) - 500)^2)))
})

test_that("the bootstrap interval is a reproducible percentile interval", {
  x <- tablet_supplier$weight_mg
  set.seed(3); a <- capability(x, lsl = 475, usl = 525, ci = "bootstrap", B = 1000)
  set.seed(3); b <- capability(x, lsl = 475, usl = 525, ci = "bootstrap", B = 1000)
  expect_equal(a$indices, b$indices)
  set.seed(3)
  xb <- matrix(sample(x, 30 * 1000, replace = TRUE), 30, 1000)
  bs <- pmin(colMeans(xb) - 475, 525 - colMeans(xb)) / (3 * apply(xb, 2, sd))
  expect_equal(a$indices$lower[a$indices$index == "Cpk"], quantile(bs, 0.025, names = FALSE))
  expect_equal(a$indices$upper[a$indices$index == "Cpk"], quantile(bs, 0.975, names = FALSE))
  expect_equal(a$lower_one_sided, quantile(bs, 0.05, names = FALSE))
})

test_that("capability validates its input", {
  x <- tablet_supplier$weight_mg
  expect_error(capability(x), "at least one specification")
  expect_error(capability(x, lsl = 525, usl = 475), "below")
  expect_error(capability(c(1, 2), lsl = 0, usl = 5), "at least 3")
  expect_error(capability(c(1, 1, 1), lsl = 0, usl = 5), "no variation")
  expect_error(capability(x, lsl = 475, usl = 525, conf = 1.5), "`conf`")
  expect_error(capability(x, lsl = 475, usl = 525, sigma = 5, ci = "bootstrap"), "sigma = NULL")
  expect_error(capability(c(1, NA, 3), lsl = 0, usl = 5), "no missing")
})

test_that("sim_capability gives the true indices and the documented sampling behavior", {
  set.seed(7)
  s <- sim_capability(30, mu = 500, sigma = 25 / (3 * 1.33), lsl = 475, usl = 525, nsim = 20000)
  expect_equal(attr(s, "true_cpk"), 1.33); expect_equal(attr(s, "true_cp"), 1.33)
  expect_equal(names(s), c("cp", "cpk", "lower", "upper", "lower_one"))
  expect_gt(mean(s$cp), 1.33)                                     # the sample Cp is biased upward (1/s)
  expect_equal(sd(s$cpk), sqrt(1 / (9 * 30) + 1.33^2 / (2 * 29)), tolerance = 0.08)   # large-sample SD
  # Cp-hat = Cp sqrt((n - 1)/chi-square): its median follows from the chi-square median
  expect_equal(median(s$cp), 1.33 * sqrt(29 / qchisq(0.5, 29)), tolerance = 0.01)
  # the interval covers the truth most of the time (about 93% to 95%)
  cover <- mean(s$lower < 1.33 & s$upper > 1.33)
  expect_gt(cover, 0.91); expect_lt(cover, 0.97)
  # one simulated sample agrees with capability()
  set.seed(8); one <- sim_capability(10, 500, 5, 475, 525, nsim = 1)
  set.seed(8); y <- rnorm(10, 500, 5)
  expect_equal(one$cpk, capability(y, lsl = 475, usl = 525)$indices$estimate[4])
  expect_equal(one$lower, capability(y, lsl = 475, usl = 525)$indices$lower[4])
})

test_that("sim_capability validates its input", {
  expect_error(sim_capability(2, 500, 5, 475, 525), "`n` must be")
  expect_error(sim_capability(30, 500, 0, 475, 525), "`sigma`")
  expect_error(sim_capability(30, 500, 5, 525, 475), "below")
  expect_error(sim_capability(30, Inf, 5, 475, 525), "`mu`")
})

test_that("dpmo and sigma_level follow the six sigma conventions", {
  expect_equal(dpmo(1), 1e6 * pnorm(-3))
  expect_equal(dpmo(2, shift = 1.5), 1e6 * pnorm(-4.5))
  expect_equal(round(dpmo(2, shift = 1.5), 1), 3.4)
  expect_equal(round(dpmo(1.33), 0), 33)
  expect_equal(sigma_level(dpmo(2, shift = 1.5)), 6)
  expect_equal(round(sigma_level(66807), 2), 3)
  expect_equal(sigma_level(dpmo(1.5, shift = 0), shift = 0), 4.5)
  expect_error(dpmo("a"), "`cpk`")
  expect_error(sigma_level(0), "between")
  expect_error(sigma_level(5, shift = Inf), "`shift`")
})

test_that("cpk_oc is the normal approximation and agrees with simulation", {
  expect_equal(cpk_oc(1.4, 30, 1.4), 0.5)
  p <- cpk_oc(c(1, 1.33, 1.67), n = 30, k = 1.4)
  expect_true(all(diff(p) > 0))
  se <- sqrt(1 / 270 + 1.33^2 / 58)
  expect_equal(p[2], pnorm((1.33 - 1.4) / se))
  set.seed(9)
  s <- sim_capability(30, 500, 25 / (3 * 1.33), 475, 525, nsim = 20000)
  expect_lt(abs(mean(s$cpk >= 1.4) - p[2]), 0.08)                  # the approximation is rough at n = 30
  expect_error(cpk_oc(-1, 30, 1.4), "positive")
})

test_that("capability_plan finds the smallest plan meeting both risks", {
  pl <- capability_plan(1.67, 1.0, alpha = 0.05, beta = 0.10)
  expect_equal(pl$n, 22); expect_equal(pl$k, 1.23, tolerance = 0.01)
  expect_equal(pl$producer_risk, 0.05, tolerance = 1e-6)
  expect_lte(pl$consumer_risk, 0.10)
  # one unit fewer does not meet the consumer's risk
  n <- pl$n - 1; k <- 1.67 - qnorm(0.95) * sqrt(1 / (9 * n) + 1.67^2 / (2 * (n - 1)))
  expect_gt(cpk_oc(1.0, n, k), 0.10)
  expect_equal(capability_plan(1.67, 1.33)$n, 95)
  expect_error(capability_plan(1.0, 1.33), "below")
  expect_error(capability_plan(1.67, 1.5, max_n = 10), "No plan")
  expect_error(capability_plan(1.67, 1.0, alpha = 2), "`alpha`")
})
