test_that("rdata returns a function", {
  model <- make_exp_model()
  rdata_fn <- rdata(model)
  expect_true(is.function(rdata_fn))
})

test_that("generated data has correct structure", {
  model <- make_exp_model()
  rdata_fn <- rdata(model)
  set.seed(42)
  df <- rdata_fn(theta = c(0.1, 0.2, 0.3), n = 100, tau = Inf, p = 0.5)

  expect_true(is.data.frame(df))
  expect_equal(nrow(df), 100)
  expect_true("t" %in% names(df))
  expect_true("omega" %in% names(df))
  expect_true("x1" %in% names(df))
  expect_true("x2" %in% names(df))
  expect_true("x3" %in% names(df))
})

test_that("all observations are exact when tau = Inf", {
  model <- make_exp_model()
  rdata_fn <- rdata(model)
  set.seed(42)
  df <- rdata_fn(theta = c(0.1, 0.2, 0.3), n = 100, tau = Inf, p = 0)

  expect_true(all(df$omega == "exact"))
})

test_that("right-censoring at tau works", {
  model <- make_exp_model()
  rdata_fn <- rdata(model)
  set.seed(42)
  df <- rdata_fn(theta = c(0.1, 0.2, 0.3), n = 500, tau = 2, p = 0)

  expect_true(any(df$omega == "right"))
  expect_true(any(df$omega == "exact"))
  # Right-censored observations have t = tau
  right_idx <- df$omega == "right"
  expect_true(all(df$t[right_idx] == 2))
  # Right-censored have empty candidate sets
  expect_true(all(!df$x1[right_idx] & !df$x2[right_idx] & !df$x3[right_idx]))
})

test_that("masking probability p controls candidate set sizes", {
  model <- make_exp_model()
  rdata_fn <- rdata(model)

  # p = 0: singleton candidate sets
  set.seed(42)
  df0 <- rdata_fn(theta = c(0.1, 0.2, 0.3), n = 500, tau = Inf, p = 0)
  cands0 <- df0$x1 + df0$x2 + df0$x3
  expect_true(all(cands0 == 1))

  # p = 1: all components in candidate set
  set.seed(42)
  df1 <- rdata_fn(theta = c(0.1, 0.2, 0.3), n = 500, tau = Inf, p = 1)
  cands1 <- df1$x1 + df1$x2 + df1$x3
  expect_true(all(cands1 == 3))
})

test_that("rdata errors on invalid parameters", {
  model <- make_exp_model()
  rdata_fn <- rdata(model)
  expect_error(rdata_fn(theta = c(-0.1, 0.2, 0.3), n = 10),
               "All parameters must be positive")
  expect_error(rdata_fn(theta = c(0.1, 0.2, 0.3), n = 10, p = -1),
               "p must be in")
  expect_error(rdata_fn(theta = c(0.1, 0.2, 0.3), n = 10, p = 2),
               "p must be in")
})

test_that("rdata works with custom column names", {
  model <- dfr_series_md(
    components = list(dfr_exponential(0.1), dfr_exponential(0.2)),
    lifetime = "time", omega = "status", candset = "c"
  )
  rdata_fn <- rdata(model)
  set.seed(42)
  df <- rdata_fn(theta = c(0.1, 0.2), n = 50, tau = Inf)
  expect_true("time" %in% names(df))
  expect_true("status" %in% names(df))
  expect_true("c1" %in% names(df))
  expect_true("c2" %in% names(df))
})

test_that("rdata for Weibull series produces valid data", {
  model <- make_weibull_model()
  rdata_fn <- rdata(model)
  set.seed(42)
  df <- rdata_fn(theta = c(2, 100, 1.5, 200), n = 50, tau = Inf, p = 0.5)
  expect_true(all(df$t > 0))
  expect_true(all(df$omega == "exact"))
})
