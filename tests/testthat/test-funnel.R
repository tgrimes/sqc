test_that("the mid-p z-score is the normal quantile of the mid-p value", {
  z <- zmid(5, 3, 1, "ratio")
  expect_equal(z, qnorm(ppois(4, 3) + 0.5 * dpois(5, 3)))
  zb <- zmid(7, 50, 0.1, "proportion")
  expect_equal(zb, qnorm(pbinom(6, 50, 0.1) + 0.5 * dbinom(7, 50, 0.1)))
  expect_lt(zmid(0, 2, 1, "ratio"), 0); expect_gt(zmid(9, 2, 1, "ratio"), 3)
})

test_that("the funnel plot flags by the z-score and its statistics are consistent", {
  f <- funnel_plot(hospital_infections$observed, hospital_infections$predicted, label = hospital_infections$hospital)
  s <- f$statistics
  expect_equal(s$value, s$count / s$size)
  expect_equal(s$z, zmid(s$count, s$size, 1, "ratio"))
  expect_equal(s$outside_998 != "", abs(s$z) > qnorm(1 - 0.001))
  expect_equal(s$outside_95 != "", abs(s$z) > qnorm(1 - 0.025))
  expect_true(all(s$outside_998[s$outside_998 != ""] == "high" | s$outside_998[s$outside_998 != ""] == "low"))
  expect_true(all(s$outside_95[s$outside_998 != ""] != ""))                    # beyond 99.8% implies beyond 95%
  expect_output(print(f), "funnel plot")
  expect_s3_class(autoplot(f), "ggplot"); expect_s3_class(as.data.frame(f), "data.frame")
})

test_that("the simulated infection data tell their story", {
  d <- hospital_infections
  f <- funnel_plot(d$observed, d$predicted, label = d$hospital)
  s <- f$statistics
  top <- which.max(s$value)
  expect_equal(s$label[top], "H028"); expect_lt(s$size[top], 1)
  expect_equal(s$outside_998[top], "")                                           # the worst SIR is inside the 99.8% limits
  expect_true(all(s$outside_998[s$label %in% c("H017", "H058", "H091", "H103")] == "high"))
  expect_equal(sum(s$outside_998 != ""), 4)
  expect_gte(sum(s$outside_95 != ""), 7)
  expect_lt(sum(d$predicted[order(-s$value)[1:10]] < 2.5) , 11)
  expect_gte(sum(d$predicted[order(-s$value)[1:10]] < 2.5), 7)                   # most of the top 10 by SIR are small hospitals
})

test_that("the funnel curves are consistent with the flags", {
  f <- funnel_plot(hospital_infections$observed, hospital_infections$predicted)
  cv <- f$curves
  up <- cv[cv$level == "99.8%", ]
  i <- which.min(abs(up$size - 8))
  ucount <- up$upper[i] * up$size[i]
  expect_gt(zmid(round(ucount), up$size[i], 1, "ratio"), qnorm(1 - 0.001))
  expect_lte(zmid(round(ucount) - 1, up$size[i], 1, "ratio"), qnorm(1 - 0.001))
  lo <- cv[cv$level == "95%", ]; j <- which.min(abs(lo$size - 10))
  expect_lt(lo$lower[j], 1); expect_gt(lo$upper[j], 1)
  w <- cv[cv$level == "99.8%", ]; k <- which.min(abs(w$size - 10))
  expect_gt(w$upper[k], lo$upper[j])                                             # the 99.8% funnel is wider
})

test_that("a funnel plot for proportions uses the binomial", {
  set.seed(2)
  n <- sample(50:800, 40); x <- rbinom(40, n, 0.12)
  f <- funnel_plot(x, n, type = "proportion")
  expect_equal(f$target, sum(x) / sum(n))
  expect_equal(f$statistics$z, zmid(x, n, sum(x) / sum(n), "proportion"))
  expect_s3_class(autoplot(f), "ggplot")
  expect_error(funnel_plot(c(5, 70), c(60, 60), type = "proportion"))
})

test_that("overdispersion over-flags, and the adjustments widen the limits", {
  set.seed(10)
  e <- 4 * hospital_infections$predicted                                           # larger hospitals (four years of data)
  m <- length(e)
  s0 <- sim_funnel(e, tau = 0, reps = 200)
  s1 <- sim_funnel(e, tau = 0.3, reps = 200)
  expect_lt(mean(s0$out_95) / m, 0.08)                                           # near the nominal 5% (discreteness)
  expect_lt(mean(s0$out_998) / m, 0.01)
  expect_gt(mean(s1$out_95) / m, 0.10)
  expect_gt(mean(s1$out_998) / m, 0.015)                                         # many more flagged than the nominal 0.2%
  s2 <- sim_funnel(e, tau = 0.3, reps = 200, overdispersion = "multiplicative")
  s3 <- sim_funnel(e, tau = 0.3, reps = 200, overdispersion = "additive")
  expect_lt(mean(s2$out_998), mean(s1$out_998) / 2); expect_lt(mean(s3$out_998), mean(s1$out_998) / 2)
  expect_lt(mean(s2$out_95) / m, 0.09); expect_lt(mean(s3$out_95) / m, 0.09)
  set.seed(11)
  th <- exp(0.4 * rnorm(m) - 0.08)
  y <- rpois(m, th * e)
  od <- funnel_plot(y, e); mu <- funnel_plot(y, e, overdispersion = "multiplicative"); ad <- funnel_plot(y, e, overdispersion = "additive")
  expect_gt(od$phi, 1.5)
  expect_gt(mu$phi, 1); expect_gt(ad$tau2, 0)
  up <- function(f, sz) { cv <- f$curves[f$curves$level == "99.8%", ]; cv$upper[which.min(abs(cv$size - sz))] }
  expect_gt(up(mu, 30), up(od, 30)); expect_gt(up(ad, 30), up(od, 30))
  expect_lte(sum(mu$statistics$outside_998 != ""), sum(od$statistics$outside_998 != ""))
  expect_lte(sum(ad$statistics$outside_998 != ""), sum(od$statistics$outside_998 != ""))
  expect_output(print(ad), "additive")
  # without extra variation the estimated factor is near 1 or below and nothing changes
  set.seed(12)
  y0 <- rpois(m, e)
  expect_equal(funnel_plot(y0, e, overdispersion = "multiplicative")$statistics$z_adj,
               funnel_plot(y0, e)$statistics$z_adj / sqrt(max(1, funnel_plot(y0, e)$phi)))
})

test_that("the arguments are checked", {
  expect_error(funnel_plot(c(1, 2), c(1, 2)))
  expect_error(funnel_plot(c(1, 2, 3), c(1, 2)))
  expect_error(funnel_plot(c(1, 2, -3), c(1, 2, 3)))
  expect_error(funnel_plot(c(1, 2, 3), c(1, 2, 0)))
  expect_error(funnel_plot(c(1, 2, 3), c(1, 2, 3), levels = c(0.99, 0.95)))
  expect_error(funnel_plot(c(1, 2, 3), c(1, 2, 3), label = c("a", "b")))
  expect_error(sim_funnel(c(1, 2, 3), type = "proportion"))
  expect_error(sim_funnel(c(1, 2, 3), tau = -1))
})

test_that("the Pareto chart sorts the categories and accumulates the percentage", {
  p <- pareto_chart(c(A = 5, B = 50, C = 20, D = 25))
  expect_s3_class(p, "ggplot")
  d <- p$data
  expect_equal(as.character(d$category), c("B", "D", "C", "A"))
  expect_equal(d$cum, c(0.5, 0.75, 0.95, 1))
  p2 <- pareto_chart(c(A = 5, B = 50, C = 20, D = 25, E = 1), top = 2)
  expect_equal(as.character(p2$data$category), c("B", "D", "Other"))
  expect_equal(p2$data$count[3], 26)
  expect_error(pareto_chart(c(1, 2, 3)))
  expect_error(pareto_chart(data.frame(a = 1)))
})
