test_that("doe_plan randomizes completely or within blocks and is reproducible", {
  p <- doe_plan(c("A", "B", "C"), replicates = 4, seed = 1)
  expect_equal(nrow(p), 12); expect_equal(p$run, 1:12); expect_equal(as.vector(table(p$treatment)), c(4, 4, 4))
  expect_equal(doe_plan(c("A", "B", "C"), replicates = 4, seed = 1), p)
  b <- doe_plan(c("A", "B"), blocks = 10, seed = 2)
  expect_equal(nrow(b), 20); expect_equal(names(b), c("run", "block", "treatment"))
  expect_true(all(tapply(b$treatment, b$block, function(x) setequal(x, c("A", "B")))))
  expect_gt(mean(b$treatment[b$run %% 2 == 1] == "A"), 0.1); expect_lt(mean(b$treatment[b$run %% 2 == 1] == "A"), 0.9)   # the first foot is not always A
  bl <- doe_plan(c("A", "B"), blocks = c("Ann", "Bob", "Cy"), replicates = 2, seed = 3)
  expect_equal(nrow(bl), 12); expect_equal(unique(bl$block), c("Ann", "Bob", "Cy"))
  set.seed(10); a <- runif(1); set.seed(10); doe_plan(c("A", "B"), seed = 5); expect_equal(runif(1), a)       # the random number state is restored
  expect_error(doe_plan("A")); expect_error(doe_plan(c("A", "A"))); expect_error(doe_plan(c("A", "B"), replicates = 0))
})

test_that("the paired test keeps its power as the units differ more and the two-sample test loses it", {
  set.seed(4)
  s <- sim_blocking(n = 10, delta = 0.4, sigma_b = c(0.01, 2), sigma_e = 0.3, reps = 3000)
  pw <- function(sb, tst) s$power[s$sigma_b == sb & s$test == tst]
  expect_equal(pw(0.01, "paired"), power.t.test(n = 10, delta = 0.4, sd = sqrt(2) * 0.3, type = "paired")$power, tolerance = 0.05)
  expect_equal(pw(2, "paired"), pw(0.01, "paired"), tolerance = 0.08)
  expect_lt(pw(2, "two-sample"), 0.12); expect_gt(pw(0.01, "two-sample"), 0.3)
  z <- sim_blocking(n = 10, delta = 0, sigma_b = 1, sigma_e = 0.3, reps = 4000)
  expect_equal(z$power, c(0.05, 0.05), tolerance = 0.25)          # false alarm probability
  expect_equal(names(s), c("sigma_b", "test", "power"))
  expect_error(sim_blocking(n = 2, delta = 1, sigma_b = 1, sigma_e = 1)); expect_error(sim_blocking(n = 10, delta = 1, sigma_b = -1, sigma_e = 1))
})

test_that("without randomization a foot effect biases the estimate and the test", {
  set.seed(6)
  no <- sim_randomization(foot_effect = 0.5, sigma_b = 2, sigma_e = 0.3, randomize = FALSE, reps = 2000)
  yes <- sim_randomization(foot_effect = 0.5, sigma_b = 2, sigma_e = 0.3, randomize = TRUE, reps = 2000)
  expect_equal(no$bias, -0.5, tolerance = 0.05)
  expect_lt(abs(yes$bias), 0.05)
  expect_gt(no$reject, 0.5); expect_lt(yes$reject, 0.12)
  expect_error(sim_randomization(foot_effect = 1, sigma_b = 1, sigma_e = 1, randomize = NA))
})

test_that("pseudoreplication inflates the false alarm probability", {
  set.seed(8)
  s <- sim_pseudorep(n_units = 5, m = 3, sigma_unit = 1, sigma_meas = 0.3, reps = 4000)
  expect_gt(s$reject[1], 0.25); expect_equal(s$reject[2], 0.05, tolerance = 0.3)
  s1 <- sim_pseudorep(n_units = 5, m = 1, sigma_unit = 1, sigma_meas = 0.3, reps = 3000)
  expect_equal(s1$reject[1], 0.05, tolerance = 0.3)                 # one measurement per unit: no pseudoreplication
  expect_error(sim_pseudorep(n_units = 1, m = 3, sigma_unit = 1, sigma_meas = 1))
})

test_that("factorial_design builds the standard order, the randomization, and fractions", {
  d <- factorial_design(c("A", "B", "C"), randomize = FALSE)
  expect_equal(nrow(d), 8); expect_equal(d$A, rep(c(-1, 1), 4)); expect_equal(d$B, rep(c(-1, -1, 1, 1), 2)); expect_equal(d$C, rep(c(-1, 1), each = 4))
  expect_equal(d$run_order, 1:8); expect_s3_class(d, "factorial_design"); expect_output(print(d), "factorial_design")
  r <- factorial_design(c("A", "B", "C"), replicates = 2, seed = 1)
  expect_equal(nrow(r), 16); expect_setequal(r$run_order, 1:16); expect_equal(as.vector(table(r$std_order)), rep(2, 8))
  expect_equal(factorial_design(c("A", "B", "C"), replicates = 2, seed = 1), r)
  expect_equal(colMeans(r[c("A", "B", "C")]), c(A = 0, B = 0, C = 0))
  l <- factorial_design(list(wing = c("short", "long"), clip = c("no", "yes")), randomize = FALSE)
  expect_equal(l$wing_level, rep(c("short", "long"), 2)); expect_equal(l$clip_level, rep(c("no", "yes"), each = 2))
  f <- factorial_design(c("A", "B", "C"), generators = list(C = c("A", "B")), randomize = FALSE)
  expect_equal(nrow(f), 4); expect_equal(f$C, f$A * f$B)
  expect_error(factorial_design("A")); expect_error(factorial_design(c("A", "A"))); expect_error(factorial_design(c("A", "B", "C"), generators = list(D = c("A", "B"))))
})

test_that("the alias structure follows the defining relation", {
  a <- alias_structure(c("A", "B", "C"), list(C = c("A", "B")))
  expect_equal(attr(a, "resolution"), 3L); expect_equal(attr(a, "defining"), "ABC")
  expect_setequal(a$aliases, c("A = B:C", "B = A:C", "C = A:B"))
  b <- alias_structure(c("A", "B", "C", "D"), list(D = c("A", "B", "C")))
  expect_equal(attr(b, "resolution"), 4L)
  expect_true("A = B:C:D" %in% b$aliases); expect_true("A:B = C:D" %in% b$aliases); expect_equal(nrow(b), 7)
  c5 <- alias_structure(c("A", "B", "C", "D", "E"), list(D = c("A", "B"), E = c("A", "C")))
  expect_equal(attr(c5, "resolution"), 3L); expect_setequal(attr(c5, "defining"), c("ABD", "ACE", "BCDE"))
  expect_error(alias_structure(c("A", "B", "C"), list()))
})

test_that("factorial_effects gives the contrasts and agrees with a linear model", {
  d <- factorial_design(c("A", "B"), replicates = 3, randomize = FALSE)
  d$y <- 10 + 2 * d$A + 1 * d$B + 3 * d$A * d$B
  fe <- factorial_effects(d, y, c("A", "B"))
  expect_equal(fe$effects$effect, c(4, 2, 6)); expect_equal(fe$intercept, 10); expect_equal(fe$effects$term, c("A", "B", "A:B"))
  fe3 <- factorial_effects(helicopter, flight_time, c("wing", "body", "clip"))
  lmf <- lm(flight_time ~ wing * body * clip, helicopter)
  expect_equal(unname(fe3$effects$effect), unname(2 * coef(lmf)[-1][c("wing", "body", "clip", "wing:body", "wing:clip", "body:clip", "wing:body:clip")]))
  expect_equal(fe3$sigma, summary(lmf)$sigma); expect_equal(fe3$df, 8)
  expect_equal(fe3$effects$p[fe3$effects$term == "wing:clip"], summary(lmf)$coefficients["wing:clip", 4])
  expect_equal(fe3$effects$se[1], 2 * fe3$sigma / sqrt(16))
  expect_equal(unname(predict(fe3, data.frame(wing = 1, body = 1, clip = 1))), unname(predict(lmf, data.frame(wing = 1, body = 1, clip = 1))))
  expect_equal(fe3$effects$term[which.min(fe3$effects$p)], "wing:clip")
  expect_output(print(fe3), "factorial_effects")
  expect_s3_class(autoplot(fe3), "ggplot"); expect_s3_class(autoplot(fe3, which = "half_normal"), "ggplot")
  # two factors with the response in a data frame with character levels
  dd <- data.frame(y = c(1, 2, 3, 5, 1.1, 2.1, 3.1, 5.2), f = rep(c("lo", "hi"), 4), g = rep(c("u", "v"), each = 2, times = 2))
  expect_equal(factorial_effects(dd, y, c("f", "g"))$effects$effect[1], mean(dd$y[dd$f == "lo"]) - mean(dd$y[dd$f == "hi"]))
})

test_that("an unreplicated experiment has no error estimate and a fraction has aliased terms", {
  d <- factorial_design(c("A", "B", "C"), randomize = FALSE)
  d$y <- 5 + 3 * d$A + 0.1 * d$B + 2 * d$A * d$C + c(0.05, -0.04, 0.02, -0.03, 0.01, 0.04, -0.02, 0)
  fe <- factorial_effects(d, y, c("A", "B", "C"))
  expect_true(is.na(fe$sigma)); expect_equal(fe$df, 0); expect_true(all(is.na(fe$effects$p)))
  expect_error(autoplot(fe)); expect_s3_class(autoplot(fe, which = "half_normal"), "ggplot")
  expect_setequal(fe$effects$term[order(-abs(fe$effects$effect))][1:2], c("A", "A:C"))
  fr <- factorial_design(c("A", "B", "C"), generators = list(C = c("A", "B")), randomize = FALSE)
  fr$y <- 5 + 3 * fr$A + 2 * fr$B + 0.5 * fr$C
  ff <- factorial_effects(fr, y, c("A", "B", "C"))
  expect_equal(ff$terms, c("A", "B", "C")); expect_setequal(ff$aliased, c("A:B", "A:C", "B:C", "A:B:C"))
  expect_error(factorial_effects(d, y, "A")); expect_error(factorial_effects(d, y, c("A", "zz")))
  expect_error(predict(fe, data.frame(A = 2, B = 1, C = 1)))
})

test_that("interaction_plot and the power functions work", {
  expect_s3_class(interaction_plot(helicopter, flight_time, "wing", "clip"), "ggplot")
  expect_error(interaction_plot(helicopter, flight_time, "wing", "zz"))
  p <- doe_power(delta = 2, sigma = 2, N = 32)
  ncp <- 2 / (2 * 2 / sqrt(32)); tc <- qt(0.975, 24)
  expect_equal(p, 1 - pt(tc, 24, ncp) + pt(-tc, 24, ncp))
  expect_gt(doe_power(2, 2, 32), doe_power(2, 2, 24)); expect_gt(doe_power(3, 2, 24), doe_power(2, 2, 24))
  tb <- doe_n(2, 2); expect_equal(tb$N, (2:12) * 8); expect_true(all(diff(tb$power) > 0)); expect_equal(tb$se, 2 * 2 / sqrt(tb$N))
  expect_equal(tb$ok, tb$power >= 0.8); expect_false(tb$ok[1]); expect_true(tb$ok[tb$replicates == 5]); expect_equal(doe_n(2, 2, power = 0.5)$ok, doe_n(2, 2)$power >= 0.5)
  expect_error(doe_power(2, 2, 8)); expect_error(doe_power(2, -1, 32))
})

test_that("one factor at a time misses the best combination when the interaction is large", {
  set.seed(11)
  s <- sim_ofat(b1 = 0.3, b2 = -0.1, b12 = 0.4, sigma = 0.01, r = 1, reps = 2000)
  expect_equal(s$p_best[s$strategy == "factorial"], 1); expect_equal(s$p_best[s$strategy == "one factor at a time"], 0)
  expect_equal(unname(attr(s, "means")), c(0.1, 0, -0.4, 0.3), tolerance = 1e-8)
  expect_equal(s$runs, c(3, 4))
  s0 <- sim_ofat(b1 = 1, b2 = 1, b12 = 0, sigma = 0.01, reps = 1000)
  expect_equal(s0$p_best, c(1, 1))
  sn <- sim_ofat(b1 = 0.3, b2 = -0.1, b12 = 0.4, sigma = 0.15, r = 2, reps = 4000)
  expect_gt(sn$p_best[2], sn$p_best[1]); expect_lt(sn$loss[2], sn$loss[1])
  expect_error(sim_ofat(1, 1, 0, sigma = 0))
})

test_that("fishbone draws a diagram", {
  g <- fishbone("Short flight time", list(Materials = c("paper", "clip"), Design = c("wing", "body", "width"), Method = c("drop height"), Measurement = c("stopwatch", "observer")))
  expect_s3_class(g, "ggplot")
  expect_error(fishbone("x", list(a = "b"))); expect_error(fishbone(c("x", "y"), list(a = "b", c = "d")))
})
