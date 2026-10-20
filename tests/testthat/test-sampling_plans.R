test_that("sampling_plan validates and prints", {
  plan <- sampling_plan(50, 1)
  expect_s3_class(plan, "sampling_plan")
  expect_output(print(plan), "n = 50")
  expect_output(print(plan), "binomial")
  expect_output(print(sampling_plan(50, 2, N = 200)), "hypergeometric")
  expect_error(sampling_plan(0, 0), "`n` must be")
  expect_error(sampling_plan(10, -1), "`c` must be")
  expect_error(sampling_plan(10, 10), "smaller than `n`")
  expect_error(sampling_plan(50, 1, N = 20), "`N` must be")
})

test_that("the binomial acceptance probability matches pbinom", {
  plan <- sampling_plan(50, 1)
  p <- c(0, 0.01, 0.05, 0.2, 1)
  expect_equal(sqc:::plan_pa(plan, p), pbinom(1, 50, p))
  expect_equal(sqc:::plan_pa(plan, 0.05), (0.95)^50 + 50 * 0.05 * 0.95^49)    # P(0) + P(1)
  expect_equal(sqc:::plan_pa(plan, 0.05), 0.2794, tolerance = 1e-3)
  expect_error(sqc:::plan_pa(plan, 1.2), "between 0 and 1")
})

test_that("the hypergeometric acceptance probability uses the lot size", {
  plan <- sampling_plan(50, 1, N = 200)
  expect_equal(sqc:::plan_pa(plan, 0.05), phyper(1, 10, 190, 50))
  # sampling a quarter of the lot reduces the variance of the defect count:
  # Pa is lower than the binomial when c is below the expected count (2.5) and higher when above
  binom <- sampling_plan(50, 1)
  expect_lt(sqc:::plan_pa(plan, 0.05), sqc:::plan_pa(binom, 0.05))
  expect_gt(sqc:::plan_pa(plan, 0.005), sqc:::plan_pa(binom, 0.005))
  # a large lot is close to the binomial
  big <- sampling_plan(50, 1, N = 1e6)
  expect_equal(sqc:::plan_pa(big, 0.05), sqc:::plan_pa(sampling_plan(50, 1), 0.05), tolerance = 1e-3)
})

test_that("oc_curve returns a tidy decreasing curve and compares plans", {
  oc <- oc_curve(sampling_plan(50, 1), p = seq(0, 0.2, by = 0.01))
  expect_s3_class(oc, "oc_curve")
  expect_equal(names(oc), c("plan", "p", "pa"))
  expect_equal(oc$pa[1], 1)
  expect_true(all(diff(oc$pa) <= 0))
  two <- oc_curve(list("A" = sampling_plan(50, 1), sampling_plan(100, 2)), p = c(0.01, 0.05))
  expect_equal(levels(two$plan), c("A", "n = 100, c = 2"))
  expect_equal(nrow(two), 4)
  expect_error(oc_curve("a"), "must be a sampling_plan")
})

test_that("autoplot builds with and without AQL and LTPD marks", {
  oc <- oc_curve(sampling_plan(50, 1))
  expect_s3_class(ggplot2::autoplot(oc), "ggplot")
  p <- ggplot2::autoplot(oc, aql = 0.01, ltpd = 0.05, title = "OC")
  expect_no_error(ggplot2::ggplot_build(p))
  expect_equal(p$labels$title, "OC")
})

test_that("plan_risks computes the producer's and consumer's risk", {
  r <- plan_risks(sampling_plan(50, 1), aql = 0.01, ltpd = 0.05)
  expect_equal(r$pa_aql, pbinom(1, 50, 0.01))
  expect_equal(r$producers_risk, 1 - pbinom(1, 50, 0.01))
  expect_equal(r$consumers_risk, pbinom(1, 50, 0.05))
  expect_error(plan_risks(sampling_plan(50, 1), aql = 0.05, ltpd = 0.01), "0 < aql < ltpd < 1")
  expect_error(plan_risks("plan", 0.01, 0.05), "created with sampling_plan")
})

test_that("AOQ, ATI, and AOQL follow the formulas", {
  plan <- sampling_plan(50, 1, N = 5000)
  p <- c(0.01, 0.05)
  pa <- sqc:::plan_pa(plan, p)
  expect_equal(aoq(plan, p), pa * p * (5000 - 50) / 5000)
  expect_equal(ati(plan, p), 50 + (1 - pa) * (5000 - 50))
  expect_equal(aoq(plan, 0), 0)
  expect_equal(ati(plan, 0), 50)                      # a lot with no defectives is always accepted
  lim <- aoql(plan)
  expect_equal(names(lim), c("p", "aoq"))
  grid <- seq(0, 0.5, by = 0.0005)
  expect_equal(unname(lim["aoq"]), max(aoq(plan, grid)))
  expect_true(all(aoq(plan, grid) <= lim["aoq"] + 1e-12))
  expect_error(aoq(sampling_plan(50, 1), 0.01), "needs the lot size")
})

test_that("find_plan returns the smallest plan that meets both risks", {
  res <- find_plan(aql = 0.01, alpha = 0.05, ltpd = 0.05, beta = 0.10)
  expect_equal(c(res$n, res$c), c(132, 3))
  expect_lte(res$producers_risk, 0.05)
  expect_lte(res$consumers_risk, 0.10)
  # brute force: no smaller sample size works for any acceptance number
  smaller <- expand.grid(n = 1:(res$n - 1), c = 0:10)
  smaller <- smaller[smaller$c < smaller$n, ]
  ok <- pbinom(smaller$c, smaller$n, 0.01) >= 0.95 & pbinom(smaller$c, smaller$n, 0.05) <= 0.10
  expect_false(any(ok))
  fin <- find_plan(aql = 0.01, ltpd = 0.05, N = 200)
  expect_lte(fin$n, 200)
  expect_lte(fin$consumers_risk, 0.10)
  expect_error(find_plan(aql = 0.05, ltpd = 0.01), "must be smaller")
  expect_error(find_plan(aql = 0.01, ltpd = 0.011, beta = 0.001, n_max = 20), "No plan")
})

test_that("sim_acceptance agrees with the exact OC curve", {
  set.seed(31)
  plan <- sampling_plan(50, 1, N = 5000)
  p <- c(0.01, 0.03, 0.05)
  sim <- sim_acceptance(plan, p, reps = 20000)
  exact <- sqc:::plan_pa(plan, p)
  expect_equal(names(sim), c("p", "pa_sim", "se"))
  expect_true(all(abs(sim$pa_sim - exact) < 4.5 * sim$se))
  bin <- sim_acceptance(sampling_plan(50, 1), 0.05, reps = 20000)
  expect_lt(abs(bin$pa_sim - pbinom(1, 50, 0.05)), 4.5 * bin$se)
  expect_error(sim_acceptance(plan, 1.5), "between 0 and 1")
})
