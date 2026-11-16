test_that("the tabular CUSUM follows its recursion", {
  x <- c(10.2, 9.7, 10.9, 11.4, 10.8, 11.6, 12.1, 9.9, 10.1)
  cs <- cusum_chart(x, target = 10, sigma = 1, k = 0.5, h = 3)
  z <- (x - 10) / 1
  up <- dn <- 0; cp <- cm <- numeric(length(x))
  for (t in seq_along(x)) { up <- max(0, up + z[t] - 0.5); dn <- max(0, dn - z[t] - 0.5); cp[t] <- up; cm[t] <- dn }
  expect_equal(cs$statistics$cplus, cp)
  expect_equal(cs$statistics$cminus, cm)
  expect_equal(cs$statistics$z, z)
  expect_equal(cs$first_signal, which(cp > 3)[1])
  expect_equal(cs$statistics$signal != "", cp > 3 | cm > 3)
})

test_that("the estimated start of the shift is the observation after the last zero", {
  x <- c(0.1, -0.2, 0.3, 1.2, 1.4, 1.1, 1.5, 1.3)
  cs <- cusum_chart(x, target = 0, sigma = 1, k = 0.5, h = 2)
  s <- cs$signals[1, ]
  zeros <- which(cs$statistics$cplus[seq_len(s$index)] == 0)
  expect_equal(s$start, max(zeros) + 1)
  expect_equal(s$direction, "up")
})

test_that("the lower sum is the mirror image of the upper sum", {
  set.seed(7)
  x <- rnorm(40, mean = c(rep(0, 20), rep(-1.2, 20)))
  a <- cusum_chart(x, target = 0, sigma = 1)
  b <- cusum_chart(-x, target = 0, sigma = 1)
  expect_equal(a$statistics$cminus, b$statistics$cplus)
  expect_equal(a$statistics$cplus, b$statistics$cminus)
  expect_equal(a$statistics$signal[40], "down")
})

test_that("a head start and a reset behave as documented", {
  x <- rep(0, 5)
  hs <- cusum_chart(x, target = 0, sigma = 1, k = 0.5, h = 4, head_start = 2)
  expect_equal(hs$statistics$cplus, c(1.5, 1, 0.5, 0, 0))
  y <- c(rep(2, 6), rep(0, 3), rep(2, 6))
  once <- cusum_chart(y, target = 0, sigma = 1, k = 0.5, h = 3)
  again <- cusum_chart(y, target = 0, sigma = 1, k = 0.5, h = 3, reset = TRUE)
  expect_equal(once$first_signal, again$first_signal)
  expect_gt(nrow(again$signals), 0)
  expect_equal(again$statistics$cplus[again$first_signal + 1], 1.5)    # restarted from 0 after the signal
})

test_that("the defaults use the data, and the arguments are checked", {
  x <- c(5, 6, 4, 7, 5, 6)
  cs <- cusum_chart(x)
  expect_equal(cs$target, mean(x)); expect_equal(cs$sigma, sd(x))
  expect_output(print(cs), "cusum_chart")
  expect_s3_class(as.data.frame(cs), "data.frame")
  expect_error(cusum_chart(c(1, NA, 3)))
  expect_error(cusum_chart(1))
  expect_error(cusum_chart(x, sigma = -1))
  expect_error(cusum_chart(x, k = 0))
  expect_error(cusum_chart(x, h = 3, head_start = 3))
  expect_s3_class(autoplot(cusum_chart(x, target = 5, sigma = 1)), "ggplot")
})

test_that("the EWMA follows its recursion and the exact limits approach the asymptotic ones", {
  x <- c(10.4, 9.8, 10.9, 11.1, 10.6, 11.3, 10.2)
  ew <- ewma_chart(x, target = 10, sigma = 0.5, lambda = 0.25, L = 2.7)
  z <- numeric(length(x)); prev <- 10
  for (t in seq_along(x)) { prev <- 0.25 * x[t] + 0.75 * prev; z[t] <- prev }
  expect_equal(ew$statistics$ewma, z)
  t <- seq_along(x)
  w <- 2.7 * 0.5 * sqrt(0.25 / 1.75 * (1 - 0.75^(2 * t)))
  expect_equal(ew$statistics$ucl, 10 + w)
  expect_equal(ew$statistics$lcl, 10 - w)
  expect_equal(ew$statistics$signal, z > 10 + w | z < 10 - w)
  as <- ewma_chart(x, target = 10, sigma = 0.5, lambda = 0.25, L = 2.7, limits = "asymptotic")
  expect_equal(as$statistics$ucl, rep(10 + 2.7 * 0.5 * sqrt(0.25 / 1.75), length(x)))
  big <- ewma_chart(rep(10, 60), target = 10, sigma = 0.5, lambda = 0.25, L = 2.7)
  expect_equal(big$statistics$ucl[60], as$statistics$ucl[1], tolerance = 1e-6)
  expect_lt(ew$statistics$ucl[1], as$statistics$ucl[1])
})

test_that("lambda = 1 gives the Shewhart chart, and the arguments are checked", {
  set.seed(3)
  x <- rnorm(30, 0, 1)
  ew <- ewma_chart(x, target = 0, sigma = 1, lambda = 1, L = 3)
  expect_equal(ew$statistics$ewma, x)
  expect_equal(ew$statistics$signal, abs(x) > 3)
  expect_error(ewma_chart(x, lambda = 0))
  expect_error(ewma_chart(x, lambda = 1.2))
  expect_error(ewma_chart(x, L = -1))
  expect_error(ewma_chart(x, limits = "wide"))
  expect_output(print(ew), "ewma_chart")
  expect_s3_class(autoplot(ew), "ggplot")
})

test_that("both charts find a 1 SD shift sooner than a Shewhart chart on average", {
  set.seed(11)
  m <- 300
  res <- replicate(m, {
    x <- c(rnorm(10), rnorm(60, mean = 1))
    c(cusum = cusum_chart(x, 0, 1, k = 0.5, h = 4.77)$first_signal,
      ewma = ewma_chart(x, 0, 1, lambda = 0.2, L = 2.859)$first_signal,
      shewhart = { s <- which(abs(x) > 3); if (length(s)) s[1] else 71L })
  })
  res[is.na(res)] <- 71L
  expect_lt(mean(res["cusum", ]), mean(res["shewhart", ]))
  expect_lt(mean(res["ewma", ]), mean(res["shewhart", ]))
})
