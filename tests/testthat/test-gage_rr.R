test_that("gage_rr reproduces the two-way aov table and the ANOVA estimators", {
  fit <- gage_rr(assay_study, value = hba1c, part = sample, operator = technician, interaction = "yes")
  tab <- summary(aov(hba1c ~ factor(sample) * technician, data = assay_study))[[1]]
  expect_equal(fit$anova$df, tab$Df)
  expect_equal(fit$anova$ss, tab$`Sum Sq`)
  expect_equal(fit$anova$ms, tab$`Mean Sq`)
  ms <- tab$`Mean Sq`
  p <- 10; o <- 3; r <- 2
  v <- fit$components$variance
  expect_equal(v[1], ms[4])                                  # repeatability
  expect_equal(v[3], (ms[2] - ms[3]) / (p * r))              # operator
  expect_equal(v[4], (ms[3] - ms[4]) / r)                    # part x operator
  expect_equal(v[6], (ms[1] - ms[3]) / (o * r))              # part
  expect_equal(v[2], max(0, v[3]) + max(0, v[4]))            # reproducibility
  expect_equal(v[5], v[1] + v[2])                            # gage R&R
  expect_equal(v[7], v[5] + v[6])
  # parts and operators are tested against the interaction; the interaction against the error
  expect_equal(fit$anova$f[1], ms[1] / ms[3])
  expect_equal(fit$anova$f[2], ms[2] / ms[3])
  expect_equal(fit$anova$f[3], ms[3] / ms[4])
  expect_equal(fit$p_interaction, tab$`Pr(>F)`[3])
})

test_that("the summary matches the definitions (%GRR, contribution, ndc, tolerance)", {
  fit <- gage_rr(assay_study, value = hba1c, part = sample, operator = technician, tolerance = 4)
  s <- fit$summary
  expect_true(fit$interaction_used)                        # p = 0.117 is below 0.25
  grr <- s$sd_grr^2; tot <- s$sd_total^2
  expect_equal(s$pct_grr, 100 * sqrt(grr / tot))
  expect_equal(s$pct_contribution, 100 * grr / tot)
  expect_equal(s$ndc, floor(sqrt(2) * s$sd_part / s$sd_grr))
  expect_equal(s$pct_tolerance, 100 * 6 * s$sd_grr / 4)
  expect_equal(s$sd_grr^2, s$sd_repeatability^2 + s$sd_reproducibility^2)
  expect_equal(sum(fit$components$percent_contribution[c(5, 6)]), 100)
  expect_equal(s$pct_grr, 17.8, tolerance = 0.01)
  expect_equal(s$ndc, 7)
  expect_output(print(fit), "%GRR")
  expect_output(print(fit), "kept in the model")
})

test_that("the interaction is pooled when its p-value is large, or when forced", {
  pooled <- gage_rr(assay_study, value = hba1c, part = sample, operator = technician, interaction = "no")
  expect_false(pooled$interaction_used)
  tab <- summary(aov(hba1c ~ factor(sample) * technician, data = assay_study))[[1]]
  pooled_ms <- (tab$`Sum Sq`[3] + tab$`Sum Sq`[4]) / (tab$Df[3] + tab$Df[4])
  expect_equal(pooled$components$variance[1], pooled_ms)
  expect_equal(pooled$components$variance[4], 0)
  expect_equal(pooled$components$variance[3], (tab$`Mean Sq`[2] - pooled_ms) / 20)
  expect_equal(pooled$anova$f[1], tab$`Mean Sq`[1] / pooled_ms)
  expect_output(print(pooled), "pooled with the error")
  # auto uses the 0.25 rule: a loose alpha pools, a strict one keeps
  expect_false(gage_rr(assay_study, value = hba1c, part = sample, operator = technician, alpha = 0.05)$interaction_used)
  expect_true(gage_rr(assay_study, value = hba1c, part = sample, operator = technician, alpha = 0.5)$interaction_used)
})

test_that("negative estimates are shown and treated as 0", {
  # three operators who agree closely: the operator variance is estimated as negative
  set.seed(4)
  d <- expand.grid(rep = 1:3, op = 1:3, part = 1:6)
  d$y <- rnorm(6, 0, 3)[d$part] + rnorm(nrow(d), 0, 0.5)
  fit <- gage_rr(d, value = y, part = part, operator = op)
  neg <- fit$negative
  expect_true(length(neg) > 0)
  idx <- match(neg[1], c("repeatability", "operator", "interaction", "part"))
  nm <- c("Repeatability", "  Operator", "  Part x operator", "Part-to-part")[idx]
  expect_lt(fit$components$variance[fit$components$component == nm], 0)
  expect_equal(fit$components$sd[fit$components$component == nm], 0)
  expect_output(print(fit), "Negative ANOVA estimate")
})

test_that("gage_rr validates its input", {
  expect_error(gage_rr(assay_study), "Supply `value`")
  expect_error(gage_rr(1:10, value = x, part = a, operator = b), "must be a data frame")
  expect_error(gage_rr(assay_study[-1, ], value = hba1c, part = sample, operator = technician), "balanced")
  expect_error(gage_rr(subset(assay_study, replicate == 1), value = hba1c, part = sample, operator = technician), "balanced")
  expect_error(gage_rr(subset(assay_study, technician == "A"), value = hba1c, part = sample, operator = technician), "2 operators")
  expect_error(gage_rr(assay_study, value = hba1c, part = sample, operator = technician, alpha = 2), "`alpha`")
  expect_error(gage_rr(assay_study, value = hba1c, part = sample, operator = technician, tolerance = -1), "`tolerance`")
  bad <- assay_study; bad$hba1c[3] <- NA
  expect_error(gage_rr(bad, value = hba1c, part = sample, operator = technician), "no missing")
})

test_that("sim_gage matches gage_rr on one simulated study", {
  set.seed(5)
  s <- sim_gage(4, 3, 2, 1, 0.5, 0.3, 0.2, nsim = 1, interaction = "yes")
  set.seed(5)
  P <- rnorm(4) * 1; O <- rnorm(3) * 0.5; PO <- matrix(rnorm(12), 4, 3) * 0.3; E <- array(rnorm(24), c(2, 4, 3)) * 0.2
  d <- expand.grid(rep = 1:2, part = 1:4, op = 1:3)
  d$y <- E[cbind(d$rep, d$part, d$op)] + P[d$part] + O[d$op] + PO[cbind(d$part, d$op)]
  fit <- gage_rr(d, value = y, part = part, operator = op, interaction = "yes")
  expect_equal(s$sigma2_e, fit$components$variance[1])
  expect_equal(s$sigma2_o, fit$components$variance[3])
  expect_equal(s$sigma2_po, fit$components$variance[4])
  expect_equal(s$sigma2_p, fit$components$variance[6])
  expect_equal(s$pct_grr, fit$summary$pct_grr)
  expect_equal(s$ndc, fit$summary$ndc)
})

test_that("sim_gage estimators are unbiased and match the F-based probability of a negative operator estimate", {
  set.seed(21)
  s <- sim_gage(10, 3, 2, sigma_p = 1.4, sigma_o = 0.2, sigma_po = 0.1, sigma_e = 0.15, nsim = 20000, interaction = "yes")
  expect_equal(mean(s$sigma2_e), 0.15^2, tolerance = 0.01)
  expect_equal(mean(s$sigma2_po), 0.1^2, tolerance = 0.05)
  expect_equal(mean(s$sigma2_o), 0.2^2, tolerance = 0.03)
  expect_equal(mean(s$sigma2_p), 1.4^2, tolerance = 0.02)
  expect_equal(names(s), c("sigma2_p", "sigma2_o", "sigma2_po", "sigma2_e", "p_interaction", "sd_grr", "pct_grr", "ndc"))
  # P(operator estimate < 0) = P(F_{o-1, (p-1)(o-1)} < 1/(1 + p r sigma_o^2/(sigma_e^2 + r sigma_po^2)))
  set.seed(22)
  s2 <- sim_gage(10, 3, 2, 1.4, 0.05, 0.1, 0.15, nsim = 30000, interaction = "yes")
  phi <- 0.05^2 / (0.15^2 + 2 * 0.1^2)
  expect_equal(mean(s2$sigma2_o < 0), pf(1 / (1 + 10 * 2 * phi), 2, 18), tolerance = 0.02)
  # the 'auto' rule changes the repeatability estimate only when the interaction is pooled
  set.seed(23)
  a <- sim_gage(10, 3, 2, 1.4, 0.2, 0, 0.15, nsim = 4000)
  expect_gt(mean(a$p_interaction > 0.25), 0.6)
})

test_that("sim_gage validates its input", {
  expect_error(sim_gage(1, 3, 2, 1, 1, 1, 1), "`p` must be")
  expect_error(sim_gage(3, 1, 2, 1, 1, 1, 1), "`o` must be")
  expect_error(sim_gage(3, 3, 1, 1, 1, 1, 1), "`r` must be")
  expect_error(sim_gage(3, 3, 2, -1, 1, 1, 1), "`sigma_p`")
  expect_error(sim_gage(3, 3, 2, 1, 1, 1, 0), "`sigma_e`")
  expect_error(sim_gage(3, 3, 2, 1, 1, 1, 1, alpha = 0), "`alpha`")
})
