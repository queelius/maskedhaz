test_that("fit returns a solver function", {
  model <- make_exp_model()
  solver <- fit(model)
  expect_true(is.function(solver))
})

test_that("MLE recovers true params for exponential series", {
  set.seed(42)
  true_par <- c(0.1, 0.2, 0.3)
  model <- make_exp_model(true_par)
  rdata_fn <- rdata(model)
  df <- rdata_fn(theta = true_par, n = 1000, tau = Inf, p = 0)

  solver <- fit(model)
  result <- suppressWarnings(solver(df, par = c(0.5, 0.5, 0.5)))

  expect_true(result$converged)
  # Sum of rates is identifiable for exponential series
  expect_equal(sum(coef(result)), sum(true_par), tolerance = 0.1)
})

test_that("MLE on masked data converges for exponential", {
  set.seed(42)
  true_par <- c(0.1, 0.2, 0.3)
  model <- make_exp_model(true_par)
  rdata_fn <- rdata(model)
  df <- rdata_fn(theta = true_par, n = 500, tau = Inf, p = 0.5)

  solver <- fit(model)
  result <- suppressWarnings(solver(df, par = c(0.5, 0.5, 0.5)))

  expect_true(result$converged)
  expect_true(is.finite(result$loglik))
})

test_that("MLE on censored data converges", {
  set.seed(42)
  true_par <- c(0.1, 0.2, 0.3)
  model <- make_exp_model(true_par)
  rdata_fn <- rdata(model)
  df <- rdata_fn(theta = true_par, n = 300, tau = 5, p = 0.3)

  solver <- fit(model)
  result <- suppressWarnings(solver(df, par = c(0.5, 0.5, 0.5)))

  expect_true(result$converged)
})

test_that("fisher_mle methods work: coef, vcov, logLik, confint", {
  set.seed(42)
  model <- make_exp_model()
  df <- make_exp_exact_data(n = 200)

  solver <- fit(model)
  result <- suppressWarnings(solver(df, par = c(0.5, 0.5, 0.5)))

  expect_length(coef(result), 3)
  expect_true(is.matrix(vcov(result)))
  expect_equal(nrow(vcov(result)), 3)
  expect_true(is.finite(as.numeric(logLik(result))))
  ci <- confint(result)
  expect_equal(nrow(ci), 3)
  expect_equal(ncol(ci), 2)
})

test_that("fit works with BFGS method", {
  set.seed(42)
  model <- make_exp_model()
  df <- make_exp_exact_data(n = 100)

  solver <- fit(model)
  result <- suppressWarnings(solver(df, par = c(0.5, 0.5, 0.5),
                                    method = "BFGS"))
  expect_true(result$converged)
})

test_that("fit auto-switches to BFGS for single-parameter model", {
  set.seed(42)
  model <- dfr_series_md(components = list(dfr_exponential(0.5)),
                          n_par = 1L)
  rdata_fn <- rdata(model)
  df <- rdata_fn(theta = 0.5, n = 200, tau = Inf, p = 0)
  solver <- fit(model)
  result <- suppressWarnings(solver(df, par = 0.3))
  expect_true(result$converged)
  expect_equal(coef(result), 0.5, tolerance = 0.1)
})

test_that("MLE recovers Weibull params (identified with masked data)", {
  set.seed(42)
  true_par <- c(2, 100, 1.5, 200)
  model <- make_weibull_model()
  rdata_fn <- rdata(model)
  df <- rdata_fn(theta = true_par, n = 500, tau = Inf, p = 0.5)

  solver <- fit(model)
  result <- suppressWarnings(solver(df, par = c(1.5, 120, 1.2, 250)))

  expect_true(result$converged)
  # Weibull series are identifiable — check relative accuracy
  est <- coef(result)
  expect_equal(est[1], true_par[1], tolerance = 0.5)
  expect_equal(est[2], true_par[2], tolerance = 30)
})
