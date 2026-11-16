test_that("Gauss-Legendre nodes integrate polynomials exactly", {
  gl <- gauss_legendre(12)
  expect_equal(sum(gl$w), 2)
  expect_equal(sum(gl$w * gl$x^2), 2 / 3)
  expect_equal(sum(gl$w * gl$x^7), 0, tolerance = 1e-12)
  expect_equal(sum(gl$w * gl$x^10), 2 / 11)
})

test_that("the Shewhart ARL is the reciprocal of the signal probability", {
  expect_equal(shewhart_arl(0, L = 3), 1 / (2 * pnorm(-3)))
  expect_equal(shewhart_arl(2, L = 3), 1 / (pnorm(-1) + pnorm(-5)))
  expect_equal(shewhart_arl(c(-1, 1)), rep(shewhart_arl(1), 2))
  expect_equal(shewhart_arl(0, L = 3), 370.4, tolerance = 1e-3)
})

test_that("the CUSUM ARL reproduces the published tabular values", {
  # Two-sided CUSUM with k = 0.5, as tabulated in Montgomery (and Hawkins and Olwell)
  expect_equal(cusum_arl(0, k = 0.5, h = 5), 465, tolerance = 2e-3)
  expect_equal(cusum_arl(0.5, k = 0.5, h = 5), 38.0, tolerance = 2e-3)
  expect_equal(cusum_arl(1, k = 0.5, h = 5), 10.4, tolerance = 5e-3)
  expect_equal(cusum_arl(0, k = 0.5, h = 4), 168, tolerance = 2e-3)
  expect_equal(cusum_arl(1, k = 0.5, h = 4), 8.38, tolerance = 2e-3)
  # one-sided chart, no head start: the same ARL after a positive shift to within the negligible lower sum
  expect_equal(cusum_arl(1, 0.5, 5, sides = 1), 10.4, tolerance = 5e-3)
  expect_equal(cusum_arl(0, 0.5, 5, sides = 1), 930.9, tolerance = 1e-3)
})

test_that("the two-sided ARL combines the two one-sided charts", {
  up <- cusum_arl(0.7, 0.5, 4, sides = 1); dn <- cusum_arl(-0.7, 0.5, 4, sides = 1)
  expect_equal(cusum_arl(0.7, 0.5, 4), 1 / (1 / up + 1 / dn))
  expect_equal(cusum_arl(-1, 0.5, 4), cusum_arl(1, 0.5, 4))      # symmetric
  expect_equal(cusum_arl(c(0, 1, 2), 0.5, 4), c(cusum_arl(0, 0.5, 4), cusum_arl(1, 0.5, 4), cusum_arl(2, 0.5, 4)))
})

test_that("a head start (fast initial response) shortens the ARL after a shift, but not much in control", {
  expect_equal(cusum_arl(1, 0.5, 5, sides = 1, head_start = 2.5), 6.35, tolerance = 5e-3)   # Lucas and Crosier
  expect_lt(cusum_arl(1, 0.5, 5, sides = 1, head_start = 2.5), cusum_arl(1, 0.5, 5, sides = 1))
  expect_gt(cusum_arl(0, 0.5, 5, sides = 1, head_start = 2.5) / cusum_arl(0, 0.5, 5, sides = 1), 0.9)
})

test_that("the EWMA ARL is correct, and lambda = 1 is a Shewhart chart", {
  expect_equal(ewma_arl(0, lambda = 0.2, L = 2.859), 370, tolerance = 2e-3)
  expect_equal(ewma_arl(c(0, 1, 2), lambda = 1, L = 3), shewhart_arl(c(0, 1, 2), L = 3), tolerance = 1e-6)
  expect_equal(ewma_arl(-1, 0.2, 2.859), ewma_arl(1, 0.2, 2.859))
  a <- ewma_arl(c(0, 0.5, 1, 2), 0.2, 2.859)
  expect_true(all(diff(a) < 0))                                   # larger shifts are found sooner
})

test_that("the limits that match an in-control ARL are found", {
  h <- find_cusum_h(370, k = 0.5)
  expect_equal(cusum_arl(0, 0.5, h), 370, tolerance = 1e-6)
  h1 <- find_cusum_h(370, k = 0.5, sides = 1)
  expect_equal(cusum_arl(0, 0.5, h1, sides = 1), 370, tolerance = 1e-6)
  expect_gt(h, h1)                                                # two sides need a wider interval for the same ARL0
  L <- find_ewma_L(370, lambda = 0.2)
  expect_equal(L, 2.859, tolerance = 1e-3)
  expect_equal(ewma_arl(0, 0.2, L), 370, tolerance = 1e-6)
})

test_that("the ARL functions check their arguments", {
  expect_error(cusum_arl(0, k = -1, h = 5))
  expect_error(cusum_arl(0, k = 0.5, h = 5, sides = 3))
  expect_error(cusum_arl("a"))
  expect_error(ewma_arl(0, lambda = 0))
  expect_error(ewma_arl(0, lambda = 1.5))
  expect_error(shewhart_arl(0, L = -3))
  expect_error(find_cusum_h(arl0 = 1))
  expect_error(find_ewma_L(arl0 = 0.5))
})

test_that("simulated run lengths agree with the integral-equation ARLs", {
  set.seed(20261021)
  s <- sim_arl("cusum", shift = 1, reps = 4000, k = 0.5, h = 5)
  expect_s3_class(s, "sqc_chart_run_lengths")
  expect_equal(arl(s), cusum_arl(1, 0.5, 5), tolerance = 0.04)
  e <- sim_arl("ewma", shift = 1, reps = 4000, lambda = 0.2, L = 2.859)
  expect_equal(arl(e), ewma_arl(1, 0.2, 2.859), tolerance = 0.04)
  w <- sim_arl("shewhart", shift = 2, reps = 4000, L = 3)
  expect_equal(arl(w), shewhart_arl(2, 3), tolerance = 0.05)
  expect_output(print(s), "ARL")
})
