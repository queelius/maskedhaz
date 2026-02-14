test_that("right-censored observations have omega = 'right'", {
  df <- make_exp_censored_data(n = 200, tau = 3)
  right_idx <- df$omega == "right"
  expect_true(any(right_idx))
  expect_true(all(df$t[right_idx] == 3))
})

test_that("loglik handles right-censored-only data", {
  model <- make_exp_model()
  # Create purely right-censored data
  df <- data.frame(
    t = c(1, 2, 3),
    omega = rep("right", 3),
    x1 = c(FALSE, FALSE, FALSE),
    x2 = c(FALSE, FALSE, FALSE),
    x3 = c(FALSE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )
  ll_fn <- loglik(model)
  # -H_sys(t) = -(0.6)*sum(t) = -0.6*6 = -3.6
  val <- ll_fn(df, par = c(0.1, 0.2, 0.3))
  expect_equal(val, -0.6 * 6, tolerance = 1e-6)
})

test_that("left-censored contributions are finite and correct sign", {
  model <- make_exp_model()
  # Create left-censored observation
  df <- data.frame(
    t = 5,
    omega = "left",
    x1 = TRUE,
    x2 = TRUE,
    x3 = FALSE,
    stringsAsFactors = FALSE
  )
  ll_fn <- loglik(model)
  val <- ll_fn(df, par = c(0.1, 0.2, 0.3))
  expect_true(is.finite(val))
  expect_true(val < 0) # log of probability is always negative
})

test_that("left-censored matches analytical for exponential", {
  model <- make_exp_model()
  par <- c(0.1, 0.2, 0.3)
  lambda_sys <- sum(par)

  # Left censored: component failed before tau, cause in {1,2}
  df <- data.frame(
    t = 5,
    omega = "left",
    x1 = TRUE,
    x2 = TRUE,
    x3 = FALSE,
    stringsAsFactors = FALSE
  )

  ll_fn <- loglik(model)
  ll_numerical <- ll_fn(df, par = par)

  # Analytical: log(lambda_c/lambda_sys * (1 - exp(-lambda_sys * tau)))
  lambda_c <- par[1] + par[2]
  ll_analytical <- log(lambda_c) + log(-expm1(-lambda_sys * 5)) - log(lambda_sys)

  expect_equal(ll_numerical, ll_analytical, tolerance = 1e-4)
})

test_that("interval-censored uses both t and t_upper", {
  model <- make_exp_model()
  par <- c(0.1, 0.2, 0.3)

  df <- data.frame(
    t = 2,
    t_upper = 5,
    omega = "interval",
    x1 = TRUE,
    x2 = FALSE,
    x3 = TRUE,
    stringsAsFactors = FALSE
  )
  ll_fn <- loglik(model)
  val <- ll_fn(df, par = par)
  expect_true(is.finite(val))
  expect_true(val < 0)
})

test_that("interval-censored matches analytical for exponential", {
  model <- make_exp_model()
  par <- c(0.1, 0.2, 0.3)
  lambda_sys <- sum(par)

  df <- data.frame(
    t = 2,
    t_upper = 5,
    omega = "interval",
    x1 = TRUE,
    x2 = TRUE,
    x3 = FALSE,
    stringsAsFactors = FALSE
  )

  ll_fn <- loglik(model)
  ll_numerical <- ll_fn(df, par = par)

  # Analytical: log(lambda_c) - lambda_sys*a + log(1-exp(-lambda_sys*(b-a))) - log(lambda_sys)
  lambda_c <- par[1] + par[2]
  a <- 2
  b <- 5
  ll_analytical <- log(lambda_c) - lambda_sys * a +
    log(-expm1(-lambda_sys * (b - a))) - log(lambda_sys)

  expect_equal(ll_numerical, ll_analytical, tolerance = 1e-4)
})

test_that("left-censored with empty candidate set is caught by validation", {
  model <- make_exp_model()
  df <- data.frame(
    t = 5, omega = "left",
    x1 = FALSE, x2 = FALSE, x3 = FALSE,
    stringsAsFactors = FALSE
  )
  ll_fn <- loglik(model)
  expect_error(ll_fn(df, par = c(0.1, 0.2, 0.3)), "non-empty candidate set")
})

test_that("interval-censored with empty candidate set is caught", {
  model <- make_exp_model()
  df <- data.frame(
    t = 2, t_upper = 5, omega = "interval",
    x1 = FALSE, x2 = FALSE, x3 = FALSE,
    stringsAsFactors = FALSE
  )
  ll_fn <- loglik(model)
  expect_error(ll_fn(df, par = c(0.1, 0.2, 0.3)), "non-empty candidate set")
})

test_that("mixed observation types work together", {
  model <- make_exp_model()
  par <- c(0.1, 0.2, 0.3)

  df <- data.frame(
    t = c(3, 5, 2, 1),
    t_upper = c(NA, NA, NA, 4),
    omega = c("exact", "right", "left", "interval"),
    x1 = c(TRUE, FALSE, TRUE, TRUE),
    x2 = c(FALSE, FALSE, TRUE, FALSE),
    x3 = c(TRUE, FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )

  ll_fn <- loglik(model)
  val <- ll_fn(df, par = par)
  expect_true(is.finite(val))
})
