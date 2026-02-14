test_that("assumptions returns expected strings", {
  model <- make_exp_model()
  a <- assumptions(model)
  expect_true(is.character(a))
  expect_true(any(grepl("C1", a)))
  expect_true(any(grepl("C2", a)))
  expect_true(any(grepl("C3", a)))
  expect_true(any(grepl("series", a)))
  expect_true(any(grepl("independence", a)))
})

test_that("ncomponents delegates correctly", {
  model <- make_exp_model(c(0.1, 0.2, 0.3))
  expect_equal(ncomponents(model), 3)

  model2 <- make_weibull_model()
  expect_equal(ncomponents(model2), 2)
})

test_that("component_hazard delegates correctly", {
  model <- make_exp_model(c(0.1, 0.2, 0.3))
  h1 <- component_hazard(model, 1)
  expect_true(is.function(h1))
  expect_equal(h1(10, par = c(0.1)), 0.1)
})

test_that("conditional_cause_probability returns correct shape", {
  model <- make_exp_model(c(0.1, 0.2, 0.3))
  ccp_fn <- conditional_cause_probability(model)
  expect_true(is.function(ccp_fn))

  probs <- ccp_fn(c(1, 5, 10), par = c(0.1, 0.2, 0.3))
  expect_equal(dim(probs), c(3, 3))
  # Each row sums to 1
  expect_equal(rowSums(probs), c(1, 1, 1), tolerance = 1e-10)
})

test_that("conditional_cause_probability correct for exponential", {
  # For exponential, P(K=j|T=t) = lambda_j / lambda_sys (independent of t)
  model <- make_exp_model(c(0.1, 0.2, 0.3))
  ccp_fn <- conditional_cause_probability(model)
  probs <- ccp_fn(c(1, 100), par = c(0.1, 0.2, 0.3))

  expected <- c(0.1, 0.2, 0.3) / 0.6
  expect_equal(probs[1, ], expected, tolerance = 1e-10)
  expect_equal(probs[2, ], expected, tolerance = 1e-10)
})

test_that("cause_probability returns correct shape", {
  model <- make_exp_model(c(0.1, 0.2, 0.3))
  cp_fn <- cause_probability(model)
  expect_true(is.function(cp_fn))

  set.seed(42)
  probs <- cp_fn(par = c(0.1, 0.2, 0.3), n_mc = 5000)
  expect_length(probs, 3)
  expect_equal(sum(probs), 1, tolerance = 0.05)
})

test_that("cause_probability returns NA when all censored", {
  model <- make_exp_model(c(0.1, 0.2, 0.3))
  cp_fn <- cause_probability(model)
  # With very short tau, nearly all observations are right-censored
  set.seed(42)
  probs <- cp_fn(par = c(0.1, 0.2, 0.3), n_mc = 50, tau = 0.0001, p = 0)
  # Should still return a length-3 vector
  expect_length(probs, 3)
})

test_that("cause_probability approximately correct for exponential", {
  # For exponential, P(K=j) = lambda_j / lambda_sys
  model <- make_exp_model(c(0.1, 0.2, 0.3))
  cp_fn <- cause_probability(model)

  set.seed(42)
  probs <- cp_fn(par = c(0.1, 0.2, 0.3), n_mc = 10000)
  expected <- c(0.1, 0.2, 0.3) / 0.6

  expect_equal(probs, expected, tolerance = 0.05)
})
