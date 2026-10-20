# Published values (rounded), as in standard SPC tables.
table_ref <- data.frame(
  n  = 2:10,
  d2 = c(1.128, 1.693, 2.059, 2.326, 2.534, 2.704, 2.847, 2.970, 3.078),
  d3 = c(0.853, 0.888, 0.880, 0.864, 0.848, 0.833, 0.820, 0.808, 0.797),
  c4 = c(0.7979, 0.8862, 0.9213, 0.9400, 0.9515, 0.9594, 0.9650, 0.9693, 0.9727),
  A2 = c(1.880, 1.023, 0.729, 0.577, 0.483, 0.419, 0.373, 0.337, 0.308),
  A3 = c(2.659, 1.954, 1.628, 1.427, 1.287, 1.182, 1.099, 1.032, 0.975),
  B3 = c(0, 0, 0, 0, 0.030, 0.118, 0.185, 0.239, 0.284),
  B4 = c(3.267, 2.568, 2.266, 2.089, 1.970, 1.882, 1.815, 1.761, 1.716),
  D3 = c(0, 0, 0, 0, 0, 0.076, 0.136, 0.184, 0.223),
  D4 = c(3.267, 2.574, 2.282, 2.114, 2.004, 1.924, 1.864, 1.816, 1.777)
)

test_that("qc_constants() reproduces the published table", {
  got <- qc_constants(2:10)
  for (col in setdiff(names(table_ref), "n")) {
    # the published values are rounded to 3 decimals
    expect_lt(max(abs(got[[col]] - table_ref[[col]])), 6e-4, label = col)
  }
})

test_that("constants have exact closed forms for n = 2", {
  expect_equal(d2(2), 2 / sqrt(pi))
  expect_equal(d3(2), sqrt(2 * (1 - 2 / pi)))
  expect_equal(c4(2), sqrt(2 / pi))
})

test_that("larger subgroup sizes agree with the tables", {
  expect_equal(d2(25), 3.931, tolerance = 1e-3)
  expect_equal(d3(25), 0.708, tolerance = 1e-3)
})

test_that("constants are vectorized and validate n", {
  expect_length(d2(2:5), 4)
  expect_s3_class(qc_constants(3:4), "tbl_df")
  expect_error(d2(1), "at least 2")
  expect_error(c4(2.5), "whole numbers")
})

test_that("moments are cached and agree with a fresh computation", {
  cache <- sqc:::.moment_cache
  if (exists("7", envir = cache)) rm("7", envir = cache)
  expect_false(exists("7", envir = cache))

  first <- sqc:::rel_range_moments(7)
  expect_true(exists("7", envir = cache))
  expect_identical(first, sqc:::compute_rel_range_moments(7))
  expect_identical(sqc:::rel_range_moments(7), first)   # second call reads the cache
})

test_that("d2, d3, and qc_constants share the same values", {
  tab <- qc_constants(c(3, 6, 12))
  expect_identical(tab$d2, d2(c(3, 6, 12)))
  expect_identical(tab$d3, d3(c(3, 6, 12)))
  m <- sqc:::rel_range_moments(6)
  expect_identical(d2(6), m[["m1"]])
  expect_identical(d3(6), sqrt(m[["m2"]] - m[["m1"]]^2))
})

test_that("repeated sizes in one call give repeated values", {
  expect_identical(d2(c(4, 4, 5, 4)), d2(c(4, 5))[c(1, 1, 2, 1)])
})
