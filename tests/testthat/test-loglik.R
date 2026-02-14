test_that("loglik returns a function", {
  model <- make_exp_model()
  ll_fn <- loglik(model)
  expect_true(is.function(ll_fn))
})

test_that("loglik is finite for valid params on exact data", {
  df <- make_exp_exact_data(n = 50, rates = c(0.1, 0.2, 0.3))
  model <- make_exp_model()
  ll_fn <- loglik(model)
  val <- ll_fn(df, par = c(0.1, 0.2, 0.3))
  expect_true(is.finite(val))
})

test_that("loglik is -Inf for non-positive params", {
  df <- make_exp_exact_data(n = 20)
  model <- make_exp_model()
  ll_fn <- loglik(model)
  expect_equal(ll_fn(df, par = c(-0.1, 0.2, 0.3)), -Inf)
  expect_equal(ll_fn(df, par = c(0, 0.2, 0.3)), -Inf)
})

test_that("loglik at true params > loglik at distant params", {
  df <- make_exp_exact_data(n = 200, rates = c(0.1, 0.2, 0.3))
  model <- make_exp_model()
  ll_fn <- loglik(model)
  ll_true <- ll_fn(df, par = c(0.1, 0.2, 0.3))
  ll_far <- ll_fn(df, par = c(1.0, 2.0, 3.0))
  expect_gt(ll_true, ll_far)
})

test_that("loglik matches analytical formula for exponential exact data", {
  df <- make_exp_exact_data(n = 100)
  par <- c(0.1, 0.2, 0.3)
  model <- make_exp_model()
  ll_fn <- loglik(model)
  ll_numerical <- ll_fn(df, par = par)
  ll_analytical <- exp_series_loglik_analytical(df, par)
  expect_equal(ll_numerical, ll_analytical, tolerance = 1e-6)
})

test_that("loglik handles mixed exact + right censored data", {
  df <- make_exp_censored_data(n = 100, tau = 5)
  model <- make_exp_model()
  ll_fn <- loglik(model)
  val <- ll_fn(df, par = c(0.1, 0.2, 0.3))
  expect_true(is.finite(val))
})

test_that("loglik matches analytical for censored exponential data", {
  df <- make_exp_censored_data(n = 100, tau = 5)
  par <- c(0.1, 0.2, 0.3)
  model <- make_exp_model()
  ll_fn <- loglik(model)
  ll_numerical <- ll_fn(df, par = par)
  ll_analytical <- exp_series_loglik_analytical(df, par)
  expect_equal(ll_numerical, ll_analytical, tolerance = 1e-6)
})

test_that("loglik works for Weibull series", {
  model <- make_weibull_model()
  rdata_fn <- rdata(model)
  set.seed(123)
  df <- rdata_fn(theta = c(2, 100, 1.5, 200), n = 50, tau = Inf)
  ll_fn <- loglik(model)
  val <- ll_fn(df, par = c(2, 100, 1.5, 200))
  expect_true(is.finite(val))
})

test_that("loglik works for mixed series", {
  model <- make_mixed_model()
  rdata_fn <- rdata(model)
  set.seed(123)
  df <- rdata_fn(theta = c(2, 100, 0.05), n = 50, tau = Inf)
  ll_fn <- loglik(model)
  val <- ll_fn(df, par = c(2, 100, 0.05))
  expect_true(is.finite(val))
})

test_that("loglik with singleton candidate sets (known cause)", {
  # Create data where candidate set = {j} for each exact obs
  set.seed(99)
  model <- make_exp_model(c(0.1, 0.2, 0.3))
  rdata_fn <- rdata(model)
  df <- rdata_fn(theta = c(0.1, 0.2, 0.3), n = 100, tau = Inf, p = 0)
  ll_fn <- loglik(model)
  val <- ll_fn(df, par = c(0.1, 0.2, 0.3))
  expect_true(is.finite(val))
})
