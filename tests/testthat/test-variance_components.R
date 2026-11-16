test_that("vc_oneway reproduces the aov table and the ANOVA estimators", {
  fit <- vc_oneway(assay_study, value = hba1c, group = sample)
  tab <- summary(aov(hba1c ~ factor(sample), data = assay_study))[[1]]
  expect_equal(fit$anova$df, tab$Df)
  expect_equal(fit$anova$ss, tab$`Sum Sq`)
  expect_equal(fit$anova$ms, tab$`Mean Sq`)
  expect_equal(fit$anova$f[1], tab$`F value`[1])
  expect_equal(fit$anova$p[1], tab$`Pr(>F)`[1])
  msb <- tab$`Mean Sq`[1]; msw <- tab$`Mean Sq`[2]
  expect_equal(fit$components$variance[2], msw)
  expect_equal(fit$components$variance[1], (msb - msw) / 6)             # balanced: n = 6 per sample
  expect_equal(fit$components$variance[3], (msb - msw) / 6 + msw)
  expect_equal(sum(fit$components$percent[1:2]), 100)
  expect_equal(fit$icc[["estimate"]], ((msb - msw) / 6) / ((msb - msw) / 6 + msw))
  expect_equal(fit$icc[["estimate"]], (fit$anova$f[1] - 1) / (fit$anova$f[1] + 5))
  expect_true(fit$balanced)
  expect_equal(fit$n0, 6)
  expect_false(fit$negative)
})

test_that("a matrix with one row per group gives the same fit as the data frame", {
  m <- matrix(assay_study$hba1c, ncol = 6, byrow = TRUE)           # the data are ordered by sample
  expect_equal(vc_oneway(m)$components, vc_oneway(assay_study, value = hba1c, group = sample)$components)
})

test_that("the confidence intervals have the documented form", {
  fit <- vc_oneway(assay_study, value = hba1c, group = sample, conf = 0.9)
  f <- fit$anova$f[1]
  fl <- f / qf(0.95, 9, 50); fu <- f * qf(0.95, 50, 9)
  expect_equal(fit$icc[["lower"]], (fl - 1) / (fl + 5))
  expect_equal(fit$icc[["upper"]], (fu - 1) / (fu + 5))
  expect_lt(fit$icc[["lower"]], fit$icc[["estimate"]])
  expect_gt(fit$icc[["upper"]], fit$icc[["estimate"]])
  ssw <- fit$anova$ss[2]
  expect_equal(fit$sd_within_ci[["lower"]], sqrt(ssw / qchisq(0.95, 50)))
  expect_equal(fit$sd_within_ci[["upper"]], sqrt(ssw / qchisq(0.05, 50)))
  # a wider level gives a wider interval
  wide <- vc_oneway(assay_study, value = hba1c, group = sample, conf = 0.99)
  expect_lt(wide$icc[["lower"]], fit$icc[["lower"]])
})

test_that("unbalanced groups use n0", {
  d <- data.frame(g = rep(c("a", "b", "c"), times = c(2, 4, 3)),
                  y = c(10, 12, 15, 14, 16, 15, 11, 9, 10))
  fit <- vc_oneway(d, value = y, group = g)
  expect_false(fit$balanced)
  N <- 9
  expect_equal(fit$n0, (N - (4 + 16 + 9) / N) / 2)
  tab <- summary(aov(y ~ g, data = d))[[1]]
  expect_equal(fit$components$variance[1], (tab$`Mean Sq`[1] - tab$`Mean Sq`[2]) / fit$n0)
})

test_that("a negative between-group estimate is flagged and reported as 0 elsewhere", {
  # group means nearly equal: MSB < MSW
  y <- rbind(c(10, 14, 12), c(11, 13, 12), c(12, 12, 12.5), c(9, 15, 12))
  fit <- vc_oneway(y)
  expect_true(fit$negative)
  expect_lt(fit$components$variance[1], 0)
  expect_equal(fit$components$sd[1], 0)
  expect_equal(fit$components$percent[1], 0)
  expect_equal(fit$icc[["estimate"]], 0)
  expect_output(print(fit), "negative")
})

test_that("vc_oneway validates its input", {
  expect_error(vc_oneway(assay_study), "supply both")
  expect_error(vc_oneway(1:10), "must be a matrix")
  expect_error(vc_oneway(matrix(1:6, 2), value = x), "only used with a data frame")
  expect_error(vc_oneway(matrix(1:6, 1)), "At least 2 groups")
  expect_error(vc_oneway(matrix(1:2, 2)), "At least one group")
  expect_error(vc_oneway(matrix(c(1, NA, 3, 4), 2)), "no missing values")
  expect_error(vc_oneway(matrix(1:6, 2), conf = 2), "`conf`")
  expect_output(print(vc_oneway(assay_study, value = hba1c, group = sample)), "ICC")
})

test_that("sim_vc estimators are unbiased and match theory", {
  set.seed(11)
  s <- sim_vc(a = 10, n = 6, sigma_a = 1.2, sigma_e = 0.25, nsim = 20000)
  expect_equal(names(s), c("msb", "msw", "f", "sigma2_a", "sigma2_e", "icc"))
  expect_equal(nrow(s), 20000)
  expect_equal(mean(s$sigma2_a), 1.44, tolerance = 0.02)
  expect_equal(mean(s$sigma2_e), 0.0625, tolerance = 0.01)
  # E(MSB) = sigma_e^2 + n sigma_a^2
  expect_equal(mean(s$msb), 0.0625 + 6 * 1.44, tolerance = 0.02)
  expect_equal(var(s$sigma2_e), 2 * 0.0625^2 / 50, tolerance = 0.05)
  # one simulated data set agrees with vc_oneway
  set.seed(5); one <- sim_vc(4, 3, 1, 1, nsim = 1)
  set.seed(5)
  a <- rnorm(4, sd = 1); y <- matrix(rnorm(12, sd = 1), 3, 4) + rep(a, each = 3)
  fit <- vc_oneway(t(y))
  expect_equal(one$msb, fit$anova$ms[1]); expect_equal(one$msw, fit$anova$ms[2])
})

test_that("the share of negative estimates equals P(F < 1/(1 + n theta))", {
  set.seed(12)
  s <- sim_vc(a = 3, n = 6, sigma_a = 0.2, sigma_e = 0.5, nsim = 40000)   # theta = 0.16
  theta <- 0.2^2 / 0.5^2
  expect_equal(mean(s$sigma2_a < 0), pf(1 / (1 + 6 * theta), 2, 15), tolerance = 0.02)
  # no group component: P(F < 1)
  set.seed(13)
  s0 <- sim_vc(a = 3, n = 6, sigma_a = 0, sigma_e = 1, nsim = 40000)
  expect_equal(mean(s0$sigma2_a < 0), pf(1, 2, 15), tolerance = 0.02)
})

test_that("sim_vc validates its input", {
  expect_error(sim_vc(1, 3, 1, 1), "`a` must be")
  expect_error(sim_vc(3, 1, 1, 1), "`n` must be")
  expect_error(sim_vc(3, 3, -1, 1), "`sigma_a`")
  expect_error(sim_vc(3, 3, 1, 0), "`sigma_e`")
  expect_error(sim_vc(3, 3, 1, 1, nsim = 0), "`nsim`")
})
