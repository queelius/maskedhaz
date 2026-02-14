test_that("score returns a function", {
  model <- make_exp_model()
  s_fn <- score(model)
  expect_true(is.function(s_fn))
})

test_that("score has correct length", {
  model <- make_exp_model()
  df <- make_exp_exact_data(n = 50)
  s_fn <- score(model)
  s <- s_fn(df, par = c(0.1, 0.2, 0.3))
  expect_length(s, 3)
})

test_that("score at MLE is near zero", {
  set.seed(42)
  model <- make_exp_model()
  df <- make_exp_exact_data(n = 500)
  solver <- fit(model)
  result <- suppressWarnings(solver(df, par = c(0.5, 0.5, 0.5)))
  s_fn <- score(model)
  s_at_mle <- s_fn(df, par = coef(result))
  expect_true(all(abs(s_at_mle) < 0.5))
})

test_that("hess_loglik returns a function", {
  model <- make_exp_model()
  h_fn <- hess_loglik(model)
  expect_true(is.function(h_fn))
})

test_that("hessian has correct dimensions", {
  model <- make_exp_model()
  df <- make_exp_exact_data(n = 50)
  h_fn <- hess_loglik(model)
  H <- h_fn(df, par = c(0.1, 0.2, 0.3))
  expect_equal(dim(H), c(3, 3))
})

test_that("hessian is negative definite at MLE", {
  set.seed(42)
  model <- make_exp_model()
  df <- make_exp_exact_data(n = 500)
  solver <- fit(model)
  result <- suppressWarnings(solver(df, par = c(0.5, 0.5, 0.5)))
  h_fn <- hess_loglik(model)
  H <- h_fn(df, par = coef(result))
  eigenvalues <- eigen(H, only.values = TRUE)$values
  expect_true(all(eigenvalues < 0))
})

test_that("score works for Weibull model", {
  model <- make_weibull_model()
  rdata_fn <- rdata(model)
  set.seed(123)
  df <- rdata_fn(theta = c(2, 100, 1.5, 200), n = 50, tau = Inf)
  s_fn <- score(model)
  s <- s_fn(df, par = c(2, 100, 1.5, 200))
  expect_length(s, 4)
  expect_true(all(is.finite(s)))
})
